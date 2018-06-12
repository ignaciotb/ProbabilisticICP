#include <icp_simple.hpp>

ICPSimple::ICPSimple(PointCloudT &cloud_ref, const Eigen::Matrix3f &tf_noise, const Eigen::Matrix3f &pcl_noise, double delta_thr){

    cloud_ref_.reset(new PointCloudT(cloud_ref));
    pcl::compute3DCentroid(*cloud_ref_, com_trg_);

    tf_noise_ = tf_noise;
    pcl_noise_ = pcl_noise;
    rms_error_ = 100;

    // Create Xi squared dist to filter out outliers in matching
    boost::math::chi_squared chi2_dist(3);
    lambda_thr_ = boost::math::quantile(chi2_dist, delta_thr);

}

void ICPSimple::getTransformMatrix(Eigen::Matrix4f& transform_matrix){
    transform_matrix = tf_mat_;
}

double ICPSimple::getRMSError(){
    return rms_error_;
}

void ICPSimple::constructKdTree(const PointCloudT::Ptr cloud_ref){
    // Construct Kd Tree with target cloud
    kdtree_.setInputCloud(cloud_ref);
}

void ICPSimple::alignStep(PointCloudT &cloud_new){

    // Find matches tuples
    std::vector<std::tuple<PointT, PointT>> matches_vec;
    if(0.16 > rms_error_){
        matches_vec = point2PlaneAssoc(cloud_new);
    }
    else{
        matches_vec = point2PointAssoc(cloud_new);
    }

    // Solve to extract tf matrix
    computeTransformationMatrix(matches_vec, cloud_new, tf_mat_);

    // Transform cloud_new based on latest tf
    pcl::transformPointCloud(cloud_new, cloud_new, tf_mat_);

    // Root Mean Square Error to measure convergence
    rms_error_ = computeRMSError(cloud_new);

    printf("RMS Error: %f \n", rms_error_);
}


std::vector<std::tuple<PointT, PointT>> ICPSimple::point2PlaneAssoc(PointCloudT &cloud_new){

    // K nearest neighbor search
    std::cout << "Point to plane association" << std::endl;

    int K = 4;
    std::vector<int> pointIdxNKNSearch(K);
    std::vector<float> pointNKNSquaredDistance(K);

    // Point-point probabilistic association
    Eigen::Vector3f error_mean;
    double pow_mhl_dist;
    PointCloudT set_AiPCL;

    // TODO: project tf_noise_ to pcl_noise subspace with the appropiate jacobian
    Eigen::Matrix3f jacobian_1 = Eigen::Matrix3f::Zero();
    jacobian_1(0,0) = tf_mat_(0,0);
    jacobian_1(1,1) = tf_mat_(0,0);
    jacobian_1(2,2) = 1;

    // Covariance matrix of error distribution: since pcls and tf have const covariances, it can be computed only once
    Eigen::Matrix3f error_cov_inv = (pcl_noise_ + jacobian_1 * tf_noise_ * jacobian_1.transpose() + pcl_noise_).inverse();

    // For every point in transformed pcl
    std::vector<std::tuple<PointT, PointT>> matches_vec;
    for(PointT point_ni: cloud_new.points){
        // Find kNN
        if(kdtree_.nearestKSearch(point_ni, K, pointIdxNKNSearch, pointNKNSquaredDistance) > 0){

            // For every potential nearest neighbor
            for(int pointTrgId: pointIdxNKNSearch){
                // Mean and cov mat of error distribution
                error_mean = Eigen::Vector3f(cloud_ref_->points[pointTrgId].x - point_ni.x,
                                             cloud_ref_->points[pointTrgId].y - point_ni.y,
                                             cloud_ref_->points[pointTrgId].z - point_ni.z);

                // Mahalanobis distance
                pow_mhl_dist = error_mean.transpose() * error_cov_inv * error_mean;

                // If Mhl dist smaller than Xi squared threshold, add to set of compatible points
                if(pow_mhl_dist < lambda_thr_){
                    set_AiPCL.points.push_back(cloud_ref_->points[pointTrgId]);
                }
            }

            // Point to Plane association: fit a plane to the points in set_Ai
            if(set_AiPCL.points.size() >= 3){    // TODO: ensure size of set_Ai >= 3

                // Weighted center of set Ai
                PointT sum_pointsj(0,0,0);
                double sum_traces = 0;
                for(PointT pointj: set_AiPCL.points){
                    sum_traces += 1/std::pow(pcl_noise_.trace(),2);
                    sum_pointsj.getArray3fMap() += pointj.getArray3fMap() * (float) 1.0/std::pow(pcl_noise_.trace(),2);
                }
                Eigen::Vector3f p_nu = Eigen::Vector3f(sum_pointsj.getArray3fMap() * 1/sum_traces);

                // Run PCA on set Ai to define normal vector to plane
                Eigen::Vector3f normal_vec = computePCAPcl(set_AiPCL);

                // Calculate distance d_p origin to plane
                double d_p = normal_vec.dot(p_nu);

                // Orthogonal projection of point_ni over the plane fitted
                Eigen::Vector3f point_ai =  Eigen::Vector3f(point_ni.getArray3fMap()) - (Eigen::Vector3f(point_ni.getArray3fMap()).dot(normal_vec) - d_p) * normal_vec;
                Eigen::Matrix3f cov_point_ai = pcl_noise_ + tf_noise_;   // TODO: add projected contributions of other noise sources

                // Error distribution between ai and ni
                // Mean and cov mat of error distribution
                error_mean = Eigen::Vector3f(point_ai(0) - point_ni.x,
                                             point_ai(1) - point_ni.y,
                                             point_ai(2) - point_ni.z);

                error_cov_inv = (cov_point_ai + jacobian_1 * tf_noise_ * jacobian_1.transpose() + pcl_noise_).inverse();

                // Mahalanobis distance
                pow_mhl_dist = error_mean.transpose() * error_cov_inv * error_mean;

                // If Mhl under threshold, create match (ai, ni)
                if(pow_mhl_dist < lambda_thr_){
                    matches_vec.push_back(std::make_tuple(PointT(point_ai(0), point_ai(1), point_ai(2)), point_ni));
                }

                set_AiPCL.points.clear();
            }
        }
    }
    printf("Size of matches vector %d \n", matches_vec.size());

    return matches_vec;
}


std::vector<std::tuple<PointT, PointT>> ICPSimple::point2PointAssoc(PointCloudT &cloud_new){

    std::cout << "Point to point association" << std::endl;

    // K nearest neighbor search
    int K = 3;
    std::vector<int> pointIdxNKNSearch(K);
    std::vector<float> pointNKNSquaredDistance(K);

    // Point-point probabilistic association
    Eigen::Vector3f error_mean;
    double pow_mhl_dist;
    std::vector<int> set_Ai;

    // TODO: project tf_noise_ to pcl_noise subspace with the appropiate jacobian
    Eigen::Matrix3f jacobian_1 = Eigen::Matrix3f::Zero();
    jacobian_1(0,0) = tf_mat_(0,0);
    jacobian_1(1,1) = tf_mat_(0,0);
    jacobian_1(2,2) = 1;

    // Covariance matrix of error distribution: since pcls and tf have const covariances, it can be computed only once
    Eigen::Matrix3f error_cov_inv = (pcl_noise_ + jacobian_1 * tf_noise_ * jacobian_1.transpose() + pcl_noise_).inverse();

    // For every point in transformed pcl
    std::vector<std::tuple<PointT, PointT>> matches_vec;
    std::vector<int>::iterator min_match_it;
    for(PointT point_ni: cloud_new.points){
        // Find kNN
        if(kdtree_.nearestKSearch(point_ni, K, pointIdxNKNSearch, pointNKNSquaredDistance) > 0){
            // For every potential nearest neighbor
            for(int pointTrgId: pointIdxNKNSearch){
                // Mean and cov mat of error distribution
                error_mean = Eigen::Vector3f(cloud_ref_->points[pointTrgId].x - point_ni.x,
                                             cloud_ref_->points[pointTrgId].y - point_ni.y,
                                             cloud_ref_->points[pointTrgId].z - point_ni.z);

                // Mahalanobis distance
                pow_mhl_dist = error_mean.transpose() * error_cov_inv * error_mean;

                // If Mhl dist smaller than Xi squared threshold, add to set of compatible points
                if(pow_mhl_dist < lambda_thr_){
                    set_Ai.push_back(pointTrgId);
                }
            }
            // The match with smallest Mhl distance is selected
            if(!set_Ai.empty()){
                min_match_it = std::min_element(std::begin(set_Ai), std::end(set_Ai));
                matches_vec.push_back(std::make_tuple(cloud_ref_->points[set_Ai.at(std::distance(std::begin(set_Ai), min_match_it))], point_ni));
                set_Ai.clear();
            }
        }
    }
    printf("Size of matches vector %d \n", matches_vec.size());

    return matches_vec;
}


Eigen::Vector3f ICPSimple::computePCAPcl(PointCloudT& set_Ai){

    // Compute the mean of the PCL
    Eigen::Vector4f com_Ai;
    pcl::compute3DCentroid(set_Ai, com_Ai);

    // Demean point cloud
    PointCloudT set_Ai_demean;
    pcl::demeanPointCloud(set_Ai, com_Ai, set_Ai_demean);

    // Compute covariance matrix
    Eigen::Matrix3f cov_mat;
    pcl::computeCovarianceMatrix(set_Ai_demean, cov_mat);

    // Extract eigenvalues and eigenvector from cov matrix
    Eigen::EigenSolver<Eigen::Matrix3f> eigenSolver;
    eigenSolver.compute(cov_mat, true);

    Eigen::MatrixXf eigenVectors = eigenSolver.eigenvectors().real();
    Eigen::VectorXf eigenvalues = eigenSolver.eigenvalues().real();

    // Return eigenvector with smallest eigenvalue
    typedef std::vector<std::pair<double, int>> PermutationIndices;
    PermutationIndices pi;
    for (unsigned int i = 0 ; i<cov_mat.cols(); i++){
        pi.push_back(std::make_pair(eigenvalues(i), i));
    }
    sort(pi.begin(), pi.end());

    return eigenVectors.col(pi[0].second).transpose();
}


void ICPSimple::computeTransformationMatrix(const std::vector<std::tuple<PointT, PointT>>& matches_vec,
                                            const PointCloudT &cloud_new,
                                            Eigen::Matrix4f& transformation_matrix){

    // Center of mass of tf point cloud
    Eigen::Vector4f com_tf;
    pcl::compute3DCentroid(cloud_new, com_tf);

    // Demean all points in the matches
    Eigen::MatrixXf trg_demean = Eigen::MatrixXf::Zero(3, matches_vec.size());
    Eigen::MatrixXf tf_demean = Eigen::MatrixXf::Zero(matches_vec.size(), 3);

    unsigned int match_cnt = 0;
    for(std::tuple<PointT, PointT> match: matches_vec){
        trg_demean.col(match_cnt) = Eigen::Vector3f(std::get<0>(match).x,
                                                    std::get<0>(match).y,
                                                    std::get<0>(match).z) - com_trg_.head(3);

        tf_demean.row(match_cnt) = (Eigen::Vector3f(std::get<1>(match).x,
                                                    std::get<1>(match).y,
                                                    std::get<1>(match).z) - com_tf.head(3)).transpose();

        match_cnt += 1;
    }

    // Assemble the correlation matrix H = source * target'
    Eigen::Matrix3f H = trg_demean * tf_demean;

    // Compute the Singular Value Decomposition
    Eigen::JacobiSVD<Eigen::Matrix3f> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix3f u = svd.matrixU();
    Eigen::Matrix3f v = svd.matrixV();

    // Compute R = U * V'
    if(u.determinant() * v.determinant() < 0) {
        for(int x = 0; x < 3; ++x) {
            v(x, 2) *= -1;
        }
    }
    Eigen::Matrix3f R = u * v.transpose();

    // Return the correct transformation
    transformation_matrix.setIdentity();

    transformation_matrix.topLeftCorner(3, 3) = R;
    const Eigen::Vector3f Rc(R * com_tf.head(3));
    transformation_matrix.block(0, 3, 3, 1) = com_trg_.head(3) - Rc;

}

double ICPSimple::computeRMSError(PointCloudT& cloud_new){

    // Find matches tuples
    std::vector<std::tuple<PointT, PointT>> matches_vec = point2PointAssoc(cloud_new);

    // Compute RMS Error
    double rmsError = 0;
    PointT diff;
    for(std::tuple<PointT, PointT> match: matches_vec){
        diff.getArray3fMap() = std::get<0>(match).getArray3fMap() - std::get<1>(match).getArray3fMap();
        rmsError += Eigen::Vector3f(diff.x, diff.y, diff.z).norm();
    }

    return std::sqrt(rmsError / matches_vec.size());
}








