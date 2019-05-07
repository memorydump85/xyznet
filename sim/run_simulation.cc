#include <cstdio>
#include <iostream>
#include <fstream>

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
        const auto &phi = mathx::to_radians(3.141569 + 4*(-16.f + scan_ix));
        std::cout << strfmt("# ray_ix:%lu scan_ix:%lu phi:%.2f\n", ray_ix_, scan_ix, phi);

        const auto &xyz = Eigen::Vector3f(cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi));
        return rot_ * xyz;
    }

    /// Return the next ray direction
    Eigen::Vector3f next() {
        const auto& dir = peek();
        ++ray_ix_;
        return dir;
    }
};

const float LidarSensorScanGenerator::START_ANGLE = mathx::to_radians(+95.f);
const float LidarSensorScanGenerator::ANGLE_RANGE = mathx::to_radians(350.f);


static
std::size_t dbg__output_marker_(
    const Eigen::Vector3f &pos,
    const Eigen::Vector3f &orientation,
    const std::size_t offset )
{
    Eigen::Matrix<float, 4, 4> arrow_head;
    arrow_head << -0.10, +0.00, +0.10, +0.00,
                  +0.10, +0.00, +0.10, -0.25,
                  -0.05, +0.00, -0.05, +0.00,
                  +1.00, +1.00, +1.00, +1.00    ;
    const auto T = mathx::xyzrph_to_matrix(pos, orientation);
    arrow_head = T * arrow_head;

    for (std::size_t i=0; i < 4; ++i) {
        std::cout << strfmt("v %.6f %.6f %.6f\n", arrow_head(0, i), arrow_head(1, i), arrow_head(2, i));
    }
    std::cout << "o marker" << offset << "\n";
    std::cout << strfmt("f %lu %lu %lu\n", offset + 0, offset + 1, offset + 3);
    std::cout << strfmt("f %lu %lu %lu\n", offset + 3, offset + 1, offset + 2);

    return 4;
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

    std::size_t voffset = 1;
    // const auto MAX_RAY_DIST = 30.f;

    // For each position in the trajectory, simulate a Lidar Scan
    for (std::size_t i=0; i < trajectory.vertices.size(); ++i) {
        const auto &pos = trajectory.vertices[i];
        const auto &orientation = trajectory.normals[i];
        voffset += dbg__output_marker_(pos, orientation, voffset);

        LidarSensorScanGenerator scan(pos, orientation);
        for (; scan.has_next(); scan.next()) {
            float dist = 2.0; //mesh.cast_ray(pos, scan.peek());
            //if (dist >= MAX_RAY_DIST) dist = 2.0;
            const Eigen::Vector3f &p = pos + dist * scan.peek();
            std::cout << strfmt("v %.6f %.6f %.6f\n", p(0), p(1), p(2));
            ++voffset;
        }
    }
}
