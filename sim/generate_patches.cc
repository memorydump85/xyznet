#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>

#include <cppkit/string.hh>
#include <cppkit/mathx.hh>
#include <cppkit/rigidbody.hh>

#include "simulator.hh"



using namespace cppkit;


static
std::size_t dbg__output_marker_(
    const Eigen::Vector3f &pos,
    const Eigen::Vector3f &orientation,
    const std::size_t offset )
{
    Eigen::Matrix<float, 4, 6> arrow_head;
    arrow_head << +0.00, -0.10, -0.10, -0.10, +0.25, +0.0,
                  +0.00, -0.10, +0.10, +0.00, +0.00, +0.0,
                  +0.00, -0.00, +0.00, +0.10, +0.00, +0.001,
                  +1.00, +1.00, +1.00, +1.00, +1.00, +1.00 ;
    const auto T = mathx::xyzrph_to_matrix(pos, orientation);
    arrow_head = T * arrow_head;

    for (int i=0; i < arrow_head.cols(); ++i) {
        std::cout << strfmt("v %.6f %.6f %.6f\n", arrow_head(0, i), arrow_head(1, i), arrow_head(2, i));
    }
    std::cout << "o marker" << offset << "\n";
    std::cout << strfmt("f %lu %lu %lu\n", offset + 0, offset + 1, offset + 4);
    std::cout << strfmt("f %lu %lu %lu\n", offset + 0, offset + 2, offset + 4);
    std::cout << strfmt("f %lu %lu %lu\n", offset + 5, offset + 3, offset + 4);

    return arrow_head.cols();
}


//-------------------------------------
class SensorNoiseModel {
//-------------------------------------
    std::random_device rand_device_;
    std::mt19937 rand_engine_;
    std::normal_distribution<float> gaussian_;

public:
    SensorNoiseModel()
        : rand_engine_(rand_device_())
        , gaussian_(0, 0.025/3)
    {}

    float sample()
        { return gaussian_(rand_engine_); }

} sensor_noise_model;


struct Pose {
    Eigen::Vector3f xyz;
    Eigen::Vector3f rph;
};


//-------------------------------------
class SensorPoseGenerator {
//-------------------------------------
    /// Rescale a uniform random sample from the range [-1, 1] to the
    /// new range [min_, max_]
    static float rescale_(const float& v, float min_, float max_)
        { return min_ + ((max_ - min_) * ((v + 1.f) / 2.f)); }

    /// Rescale a uniform random sample from the range [-1, 1] to the
    /// new range [min_, max_]
    static Eigen::Vector3f rescale_(const Eigen::Vector3f& v, float min_, float max_)
        { return min_ + ((max_ - min_) * ((v.array() + 1.f) / 2.f)); }

    static Eigen::Vector3f polar_to_xyz_(const Eigen::Vector2f &p)
        { return Eigen::Vector3f(p(0)*cos(p(1)), p(0)*sin(p(1)), 0.5); }

public:
    std::pair<Pose, Pose> next_pose_pair() const {
        Eigen::Vector2f q = Eigen::Vector2f::Random();
        q(0) = rescale_(q(0), 3, 17);       // between 3 and 17 meters
        q(1) = rescale_(q(1), -M_PI, M_PI); // -180 to 180 degrees

        Eigen::Vector3f rph = Eigen::Vector3f::Random();
        using cppkit::mathx::to_radians;
        rph(0) = 0;                             // no roll
        rph(1) = 0;                             // no pitch
        rph(2) = rescale_(rph(2), -M_PI, M_PI); // any yaw in [-180, 180] degrees

        // Move a little, but don't change angle.
        Eigen::Vector2f q_delta = Eigen::Vector2f::Random();
        q_delta(1) = 0;

        // rph +/- 5 degrees
        Eigen::Vector3f rph_delta = Eigen::Vector3f::Random();
        rph_delta = rescale_(rph_delta.array(), to_radians(-5.f), to_radians(+5.f));

        return std::make_pair(
                    Pose{ polar_to_xyz_(q),             rph             },
                    Pose{ polar_to_xyz_(q + q_delta),   rph + rph_delta }
                );
    }

} sensor_pose_generator;


int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "  USAGE: %s <model_file.obj>\n", argv[0]);
        return -1;
    }

    std::ifstream infile(argv[1]);
    CHECK_MSG( infile.is_open(),
        cppkit::strfmt("Cannot open file: \"%s\"", argv[1]) );

    IndexedTriangleMesh mesh = IndexedTriangleMesh::from_wavefront_obj(infile, 125);
    printf("# Loaded model with %lu vertices, %lu faces\n",
        mesh.vertices.size(), mesh.vix_tuples.size());

    std::size_t vx_offset = 1;
    const auto MAX_RAY_DIST = 100.f;

    // Simulate LIDAR scans for a bunch of random pose pairs
    for (std::size_t i=0; i < 1024; ++i) {
        const auto &[pose1, pose2] = sensor_pose_generator.next_pose_pair();

        for (const auto &pose : {pose1, pose2}) {
            std::cout << strfmt("\n\n# trajectory position #%lu\n", i);
            vx_offset += dbg__output_marker_(pose.xyz, pose.rph, vx_offset);
            std::cout << '\n';

            for (LidarSensorScanGenerator scan(pose.xyz, pose.rph); scan.has_next(); scan.next()) {
                const float dist = mesh.cast_ray(pose.xyz, scan.peek()) + sensor_noise_model.sample();
                if (dist >= MAX_RAY_DIST) continue;
                const Eigen::Vector3f &p = pose.xyz + dist * scan.peek();
                std::cout << strfmt("v %.6f %.6f %.6f\n", p(0), p(1), p(2));
                ++vx_offset;
            }
        }
    }
}
