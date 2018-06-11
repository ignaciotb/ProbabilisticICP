#include <icp_simple.hpp>

ICPSimple::ICPSimple(PointCloudT &cloud_trg, const Eigen::Matrix3d& tf_noise, const Eigen::Matrix3d& pcl_noise, double delta_thr){

    cloud_trg_.reset(new PointCloudT(cloud_trg));
    pcl::compute3DCentroid(*cloud_trg_, com_trg_);

    tf_noise_ = tf_noise;
    pcl_noise_ = pcl_noise;

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

void ICPSimple::constructKdTree(const PointCloudT::Ptr cloud_trg){
    // Construct Kd Tree with target cloud
    kdtree_.setInputCloud(cloud_trg);
}

void ICPSimple::alignStep(PointCloudT &cloud_tf){

    // Find matches tuples
    std::vector<std::tuple<PointT, PointT>> matches_vec = point2PointAssoc(cloud_tf);

    // Solve to extract tf matrix
    computeTransformationMatrix(matches_vec, cloud_tf, tf_mat_);

    // Transform cloud_tf based on latest tf
    pcl::transformPointCloud(cloud_tf, cloud_tf, tf_mat_);

    // Root Mean Square Error to measure convergence
    rms_error_ = computeRMSError(cloud_tf);

    printf("RMS Error: %f \n", rms_error_);
}


std::vector<std::tuple<PointT, PointT>> ICPSimple::point2PlaneAssoc(PointCloudT &cloud_tf){

    // K nearest neighbor search
    int K = 3;
    std::vector<int> pointIdxNKNSearch(K);
    std::vector<float> pointNKNSquaredDistance(K);

    // Point-point probabilistic association
    Eigen::Vector3d error_mean;
    double pow_mhl_dist;
    PointCloudT set_AiPCL;

    // TODO: project tf_noise_ to pcl_noise subspace with the appropiate jacobian
    Eigen::Matrix3d jacobian_1 = Eigen::Matrix3d::Zero();
    jacobian_1(0,0) = tf_mat_(0,0);
    jacobian_1(1,1) = tf_mat_(0,0);
    jacobian_1(2,2) = 1;

    // Covariance matrix of error distribution: since pcls and tf have const covariances, it can be computed only once
    Eigen::Matrix3d error_cov_inv = (pcl_noise_ + jacobian_1 * tf_noise_ * jacobian_1.transpose() + pcl_noise_).inverse();

    // For every point in transformed pcl
    std::vector<std::tuple<PointT, PointT>> matches_vec;
    for(PointT point_i: cloud_tf.points){
        // Find kNN
        if(kdtree_.nearestKSearch(point_i, K, pointIdxNKNSearch, pointNKNSquaredDistance) > 0){
            // For every potential nearest neighbor
            for(int pointTrgId: pointIdxNKNSearch){
                // Mean and cov mat of error distribution
                error_mean = Eigen::Vector3d(cloud_trg_->points[pointTrgId].x - point_i.x,
                                             cloud_trg_->points[pointTrgId].y - point_i.y,
                                             cloud_trg_->points[pointTrgId].z - point_i.z);

                // Mahalanobis distance
                pow_mhl_dist = error_mean.transpose() * error_cov_inv * error_mean;

                // If Mhl dist smaller than Xi squared threshold, add to set of compatible points
                if(pow_mhl_dist < lambda_thr_){
                    set_AiPCL.points.push_back(cloud_trg_->points[pointTrgId]);
                }
            }
            // Point to Plane association: fit a plane to the points in set_Ai
            if(!set_AiPCL.points.size() >= 3){    // TODO: ensure size of set_Ai >= 3

                // Weighted center of set of points Ai
                PointT sum_pointsj(0,0,0);
                double sum_traces = 0;
                for(PointT pointj: set_AiPCL.points){
                    sum_traces += 1/std::pow(pcl_noise_.trace(),2);
                    sum_pointsj.getArray3fMap() += pointj.getArray3fMap() * (float) 1.0/std::pow(pcl_noise_.trace(),2);
                }
                Eigen::Vector3f p_nu = Eigen::Vector3f(sum_pointsj.getArray3fMap() * 1/sum_traces);

                // Normal (unitary) vector of the plane from PCA analysis





                set_AiPCL.points.clear();
            }
        }
    }
    printf("Size of matches vector %d \n", matches_vec.size());

    return matches_vec;
}


void ICPSimple::computePCAPcl(PointCloudT& set_Ai){

    // Compute the mean of the PCL
    Eigen::Vector4f com_Ai;
    pcl::compute3DCentroid(set_Ai, com_Ai);

    // Demean point cloud
    PointCloudT set_Ai_demean;
    pcl::demeanPointCloud(set_Ai, com_Ai, set_Ai_demean);

    // Compute covariance matrix
    Eigen::Matrix3d cov_mat;
    pcl::computeCovarianceMatrix(set_Ai_demean, cov_mat);

    // Extract eigenvalues and eigenvector from cov matrix
    Eigen::EigenSolver<Eigen::Matrix3d> eigenSolver;
    eigenSolver.compute(cov_mat, true);
    auto eigenvalues = eigenSolver.eigenvalues();
    auto eigenVectors = eigenSolver.eigenvectors();
}


std::vector<std::tuple<PointT, PointT>> ICPSimple::point2PointAssoc(PointCloudT &cloud_tf){

    // K nearest neighbor search
    int K = 3;
    std::vector<int> pointIdxNKNSearch(K);
    std::vector<float> pointNKNSquaredDistance(K);

    // Point-point probabilistic association
    Eigen::Vector3d error_mean;
    double pow_mhl_dist;
    std::vector<int> set_Ai;

    // TODO: project tf_noise_ to pcl_noise subspace with the appropiate jacobian
    Eigen::Matrix3d jacobian_1 = Eigen::Matrix3d::Zero();
    jacobian_1(0,0) = tf_mat_(0,0);
    jacobian_1(1,1) = tf_mat_(0,0);
    jacobian_1(2,2) = 1;

    // Covariance matrix of error distribution: since pcls and tf have const covariances, it can be computed only once
    Eigen::Matrix3d error_cov_inv = (pcl_noise_ + jacobian_1 * tf_noise_ * jacobian_1.transpose() + pcl_noise_).inverse();

    // For every point in transformed pcl
    std::vector<std::tuple<PointT, PointT>> matches_vec;
    std::vector<int>::iterator min_match_it;
    for(PointT point_i: cloud_tf.points){
        // Find kNN
        if(kdtree_.nearestKSearch(point_i, K, pointIdxNKNSearch, pointNKNSquaredDistance) > 0){
            // For every potential nearest neighbor
            for(int pointTrgId: pointIdxNKNSearch){
                // Mean and cov mat of error distribution
                error_mean = Eigen::Vector3d(cloud_trg_->points[pointTrgId].x - point_i.x,
                                             cloud_trg_->points[pointTrgId].y - point_i.y,
                                             cloud_trg_->points[pointTrgId].z - point_i.z);

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
                matches_vec.push_back(std::make_tuple(cloud_trg_->points[set_Ai.at(std::distance(std::begin(set_Ai), min_match_it))], point_i));
                set_Ai.clear();
            }
        }
    }
    printf("Size of matches vector %d \n", matches_vec.size());

    return matches_vec;
}


void ICPSimple::computeTransformationMatrix(const std::vector<std::tuple<PointT, PointT>>& matches_vec,
                                            const PointCloudT &cloud_tf,
                                            Eigen::Matrix4f& transformation_matrix){

    // Center of mass of tf point cloud
    Eigen::Vector4f com_tf;
    pcl::compute3DCentroid(cloud_tf, com_tf);

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

double ICPSimple::computeRMSError(PointCloudT& cloud_tf){

    // Find matches tuples
    std::vector<std::tuple<PointT, PointT>> matches_vec = point2PointAssoc(cloud_tf);

    // Compute RMS Error
    double rmsError = 0;
    PointT diff;
    for(std::tuple<PointT, PointT> match: matches_vec){
        diff.getArray3fMap() = std::get<0>(match).getArray3fMap() - std::get<1>(match).getArray3fMap();
        rmsError += Eigen::Vector3f(diff.x, diff.y, diff.z).norm();
    }

    return std::sqrt(rmsError / matches_vec.size());
}








