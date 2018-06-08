#ifndef ICP_SIMPLE_CPP
#define ICP_SIMPLE_CPP

#include <pcl/point_types.h>
#include <pcl/registration/icp.h>

#include <tuple>
#include <algorithm>
#include <functional>

#include <eigen3/Eigen/Core>

typedef pcl::PointXYZ PointT;
typedef pcl::PointCloud<PointT> PointCloudT;

class ICPSimple{

public:


    ICPSimple(PointCloudT &cloud_trg);

    void constructKdTree(const PointCloudT::Ptr cloud_trg);

    void alignStep(PointCloudT &cloud_tf);

    void getTransformMatrix(Eigen::Matrix4f& transform_matrix);

    double getRMSError();

private:

    // Inputs
    PointCloudT::Ptr cloud_trg_;
    Eigen::Matrix4f tf_mat_;

    // KdTree of target cloud
    pcl::KdTreeFLANN<pcl::PointXYZ> kdtree_;

    // Center of mass of target cloud
    Eigen::Vector4f com_trg_;

    // Convergence error
    double rms_error_;

    // Methods
    void computeTransformationMatrix(const std::vector<std::tuple<PointT, PointT>>& matches_vec,
                                     const PointCloudT &cloud_tf,
                                     Eigen::Matrix4f& transformation_matrix);

    std::vector<std::tuple<PointT, PointT>> matchPointClouds(PointCloudT &cloud_tf);

    double computeRMSError(PointCloudT &cloud_tf);
};

#endif // ICP_SIMPLE_CPP
