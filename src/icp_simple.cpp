#include <icp_simple.hpp>

ICPSimple::ICPSimple(PointCloudT &cloud_trg, const Eigen::Matrix4f& init_tf_mat){

    *cloud_trg_ = cloud_trg;
    tf_mat_ = init_tf_mat;
    pcl::compute3DCentroid(*cloud_trg_, com_trg_);

}

void ICPSimple::constructKdTree(const PointCloudT::Ptr cloud_trg){
    // Construct Kd Tree with target cloud
    kdtree_.setInputCloud(cloud_trg);
}

void ICPSimple::alignStep(PointCloudT &cloud_tf){

    // Transform cloud_tf based on latest tf
    pcl::transformPointCloud(cloud_tf, cloud_tf, tf_mat_);
    std::cout << "PCL transformed " << std::endl;

    // Find matches tuples
    std::vector<std::tuple<PointT, PointT>> matches_vec = matchPointClouds(cloud_tf);
    std::cout << "Correspondences found " << std::endl;

    // Solve to extract tf matrix
    computeTransformationMatrix(matches_vec, cloud_tf, tf_mat_);
    std::cout << "Transformation extracted " << std::endl;
}

void ICPSimple::getTransformMatrix(Eigen::Matrix4f& transform_matrix){
    transform_matrix = tf_mat_;
}

std::vector<std::tuple<PointT, PointT>> ICPSimple::matchPointClouds(PointCloudT &cloud_tf){

    // K nearest neighbor search
    int K = 1;
    double dmax_squared = std::pow(0.1, 2);

    std::vector<int> pointIdxNKNSearch(K);
    std::vector<float> pointNKNSquaredDistance(K);

    std::vector<std::tuple<PointT, PointT>> matches_vec;
    for(PointT point_i: cloud_tf.points){
//        count = kdtree.radiusSearch (searchPoint, radius, pointIdxRadiusSearch, pointsSquaredDistRadius); // TODO: Try radius search?
        if(kdtree_.nearestKSearch (point_i, K, pointIdxNKNSearch, pointNKNSquaredDistance) > 0){
            if(pointNKNSquaredDistance[0] <= dmax_squared){
                matches_vec.push_back(std::make_tuple(cloud_trg_->points[pointIdxNKNSearch[0]], point_i));
            }
        }
    }

    return matches_vec;
}


void ICPSimple::computeTransformationMatrix(const std::vector<std::tuple<PointT, PointT>>& matches_vec,
                                            const PointCloudT &cloud_tf,
                                            Eigen::Matrix4f& transformation_matrix){

    // Center of mass of tf point cloud
    Eigen::Vector4f com_tf;
    pcl::compute3DCentroid(cloud_tf, com_tf);

    // Demean all points in the matches
    Eigen::MatrixXf tf_demean = Eigen::MatrixXf(4, matches_vec.size());
    Eigen::MatrixXf trg_demean = Eigen::MatrixXf(matches_vec.size(), 4);

    unsigned int match_cnt = 0;
    for(std::tuple<PointT, PointT> match: matches_vec){
        trg_demean.col(match_cnt) = Eigen::Vector4f(std::get<0>(match).x,
                                                    std::get<0>(match).y,
                                                    std::get<0>(match).z,
                                                    1) - com_trg_;

        tf_demean.row(match_cnt) = Eigen::Vector4f(std::get<1>(match).x,
                                                    std::get<1>(match).y,
                                                    std::get<1>(match).z,
                                                    1) - com_tf;
        match_cnt += 1;
    }

    // Assemble the correlation matrix H = source * target'
    Eigen::Matrix3f H = (trg_demean * tf_demean).topLeftCorner(3, 3);

    // Compute the Singular Value Decomposition
    Eigen::JacobiSVD<Eigen::Matrix3f > svd (H, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix3f u = svd.matrixU();
    Eigen::Matrix3f v = svd.matrixV();

    // Compute R = V * U'
    if(u.determinant() * v.determinant() < 0) {
        for(int x = 0; x < 3; ++x) {
            v(x, 2) *= -1;
        }
    }

    Eigen::Matrix3f R = v * u.transpose();

    // Return the correct transformation
    transformation_matrix.setIdentity();

    transformation_matrix.topLeftCorner(3, 3) = R;
    const Eigen::Vector3f Rc(R * com_tf.head(3));
    transformation_matrix.block(0, 3, 3, 1) = com_trg_.head(3) - Rc;

}
