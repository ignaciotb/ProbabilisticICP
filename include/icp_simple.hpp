#ifndef ICP_SIMPLE_CPP
#define ICP_SIMPLE_CPP

#include <pcl/point_types.h>
#include <pcl/registration/icp.h>

#include <tuple>
#include <algorithm>
#include <functional>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>

#include <boost/math/distributions.hpp>

typedef pcl::PointXYZ PointT;
typedef pcl::PointCloud<PointT> PointCloudT;

class ICPSimple{

public:


    ICPSimple(PointCloudT &cloud_trg, const Eigen::Matrix3f& tf_noise, const Eigen::Matrix3f& pcl_noise, double delta_thr);

    void constructKdTree(const PointCloudT::Ptr cloud_trg);

    void alignStep(PointCloudT &cloud_tf);

    void getTransformMatrix(Eigen::Matrix4f& transform_matrix);

    double getRMSError();

    bool converged();

private:

    // Inputs
    PointCloudT::Ptr cloud_ref_;
    Eigen::Matrix3f tf_noise_;
    Eigen::Matrix3f pcl_noise_;

    // Estimated tf
    Eigen::Matrix4f tf_mat_;

    // KdTree of target cloud
    pcl::KdTreeFLANN<pcl::PointXYZ> kdtree_;

    // Center of mass of target cloud
    Eigen::Vector4f com_trg_;

    // Convergence error
    double rms_error_;
    double rms_error_prev_;

    // Aux
    double lambda_thr_;

    // Methods
    void computeTransformationMatrix(const std::vector<std::tuple<PointT, PointT>>& matches_vec,
                                     const PointCloudT &cloud_tf,
                                     Eigen::Matrix4f& transformation_matrix);

    std::vector<std::tuple<PointT, PointT>> point2PointAssoc(PointCloudT &cloud_tf);

    std::vector<std::tuple<PointT, PointT>> point2PlaneAssoc(PointCloudT &cloud_tf);

    double computeRMSError(PointCloudT &cloud_tf);

    Eigen::Vector3f computePCAPcl(PointCloudT &set_Ai);

};

#endif // ICP_SIMPLE_CPP
