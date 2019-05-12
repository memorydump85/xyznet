#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>

#include <cppkit/mathx.hh>
#include <cppkit/rigidbody.hh>

#include "simulator.hh"



using namespace cppkit;

//-------------------------------------
struct LidarSensorScanGenerator {
//-------------------------------------
    static const float START_ANGLE;
    static const float ANGLE_RANGE;
    static const std::size_t NUM_RAYS_PER_SCAN = 1080;
    static const std::size_t NUM_SCANS = 32;

    const Eigen::Vector3f pos;
    const Eigen::Matrix3f rot_;
    std::size_t ray_ix_ = 0;

    LidarSensorScanGenerator(
        const Eigen::Vector3f &xyz,
        const Eigen::Vector3f &rph )
        : pos(xyz)
        , rot_(mathx::rotation_matrix(rph(0), rph(1), rph(2)).block(0, 0, 3, 3))
    {}

    bool has_next() const {
        return ray_ix_ < (NUM_SCANS * NUM_RAYS_PER_SCAN);
    }

    /// Peek the next ray direction
    Eigen::Vector3f peek() const {
        const auto &DEL_THETA = ANGLE_RANGE / NUM_RAYS_PER_SCAN;
        const auto &theta = START_ANGLE + (ray_ix_ % NUM_RAYS_PER_SCAN) * DEL_THETA;

        const auto &scan_ix = ray_ix_ / NUM_RAYS_PER_SCAN;
        const auto &phi = mathx::to_radians(90 + 3*(-16.f + scan_ix));

        const auto &xyz = Eigen::Vector3f(cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi));
        return rot_ * xyz;
    }

    /// Return the next ray direction
    Eigen::Vector3f next() {
        const auto& dir = peek();
        ++ray_ix_;
        return dir;
    }

    /// Get a vector containing rays for an entire scan
    static
    std::vector<Eigen::Vector3f> get_scan(
        const Eigen::Vector3f &pos,
        const Eigen::Vector3f &orientation )
    {
        std::vector<Eigen::Vector3f> rays;
        rays.reserve(LidarSensorScanGenerator::NUM_SCANS * LidarSensorScanGenerator::NUM_RAYS_PER_SCAN);

        for (LidarSensorScanGenerator s(pos, orientation); s.has_next(); s.next()) {
            rays.push_back(s.peek());
        }
        return std::move(rays);
    }
};

const float LidarSensorScanGenerator::START_ANGLE = mathx::to_radians(+185.f);
const float LidarSensorScanGenerator::ANGLE_RANGE = mathx::to_radians(350.f);


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


int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "  USAGE: %s <file.obj> <trajectory.obj>\n", argv[0]);
        return -1;
    }

    const char *model_filename = argv[1];
    const char *trajectory_filename = argv[2];

    std::ifstream infile1(model_filename);
    TriangleMesh mesh = TriangleMesh::from_wavefront_obj(infile1);
    printf("# Loaded model with %lu vertices, %lu faces\n", mesh.vertices.size(), mesh.indices.size() / 3);

    // The trajectory file is also an OBJ file. We treat the vertices as
    // sensor position specifications and use the corresponding vertex
    // normals as the (roll, pitch, yaw) vector specifying the sensor
    // orientation.
    std::ifstream infile2(trajectory_filename);
    TriangleMesh trajectory = TriangleMesh::from_wavefront_obj(infile2);
    printf("# Simulating trajectory with %lu positions\n", trajectory.vertices.size());

    CHECK_MSG( trajectory.vertices.size() == trajectory.normals.size(),
                "Trajectory does not specify an orientation for all positions" );

    std::size_t vx_offset = 1;
    const auto MAX_RAY_DIST = 30.f;

    /// Noise model
    std::random_device rand_device;
    std::mt19937 rand_engine(rand_device());
    std::normal_distribution<double> noise_model_(0, 0.025/3);

    // For each position in the trajectory, simulate a Lidar Scan
    for (std::size_t i=0; i < trajectory.vertices.size(); ++i) {
        const auto &pos = trajectory.vertices[i];
        const auto &orientation = trajectory.normals[i];

        std::cout << strfmt("\n\n# trajectory position #%lu\n", i);
        vx_offset += dbg__output_marker_(pos, orientation, vx_offset);
        std::cout << '\n';

        const auto &scan_rays = LidarSensorScanGenerator::get_scan(pos, orientation);

        for (std::size_t i=0; i < scan_rays.size(); ++i) {
            const float dist = mesh.cast_ray(pos, scan_rays[i]) + noise_model_(rand_engine);
            if (dist >= MAX_RAY_DIST) continue;
            const Eigen::Vector3f &p = pos + dist * scan_rays[i];
            std::cout << strfmt("v %.6f %.6f %.6f\n", p(0), p(1), p(2));
            ++vx_offset;
        }
    }
}
