#include <cstdio>
#include <fstream>

#include <cppkit/mathx.hh>
#include <cppkit/rigidbody.hh>
#include <cppkit/test_framework.hh>

#include "simulator.hh"



using Vec3 = Eigen::Vector3f;


//
// Ray triangle
//


const Vec3 TRIANGLE[3] = {
    Vec3(+0.0f, +0.0f, +0.0f),
    Vec3(+1.0f, +0.0f, +0.0f),
    Vec3(+0.0f, +1.0f, +0.0f),
};


float intersect_ray_with_triangle(
    const Eigen::Vector3f &origin,  /// ray origin
    const Eigen::Vector3f &dir,     /// ray direction
    const Eigen::Vector3f &v0,      /// triangle vertex1
    const Eigen::Vector3f &v1,      /// triangle vertex2
    const Eigen::Vector3f &v2       /// triangle vertex3
);


DEFINE_TEST(ray_triangle_simple_intersection) {
    auto dist = intersect_ray_with_triangle(
        Vec3(+0.5f, +0.5f, +1.0f),
        Vec3(+0.0f, +0.0f, -1.0f),
        TRIANGLE[0], TRIANGLE[1], TRIANGLE[2]
    );
    CHECK(cppkit::mathx::is_close(dist, 1.0f));
}


DEFINE_TEST(ray_triangle_no_intersection) {
    // Wrong ray direction
    auto dist = intersect_ray_with_triangle(
        Vec3(+0.5f, +0.5f, -1.0f),
        Vec3(+0.0f, +0.0f, -1.0f),
        TRIANGLE[0], TRIANGLE[1], TRIANGLE[2]
    );
    CHECK(false == std::isfinite(dist));

    // Wrong ray direction
    dist = intersect_ray_with_triangle(
        Vec3(+0.5f, +0.5f, +1.0f),
        Vec3(+0.0f, +0.0f, +1.0f),
        TRIANGLE[0], TRIANGLE[1], TRIANGLE[2]
    );
    CHECK(false == std::isfinite(dist));

    // Ray intersects with triangle plane, but outside triangle
    dist = intersect_ray_with_triangle(
        Vec3(+1.0f, +1.0f, +1.0f),
        Vec3(+0.0f, +0.0f, -1.0f),
        TRIANGLE[0], TRIANGLE[1], TRIANGLE[2]
    );
    CHECK(false == std::isfinite(dist));

    // Ray is parallel to triangle
    dist = intersect_ray_with_triangle(
        Vec3(+0.5f, +0.5f, +1.0f),
        Vec3(+1.0f, +1.0f, +0.0f),
        TRIANGLE[0], TRIANGLE[1], TRIANGLE[2]
    );
    CHECK(false == std::isfinite(dist));
}


//
// Ray box
//


const Vec3 BOX[2] = {
    Vec3(-1.0f, -1.0f, -1.0f),
    Vec3(+1.0f, +1.0f, +1.0f),
};


bool ray_intersects_box(
    const Eigen::Vector3f &origin,  /// ray origin
    const Eigen::Vector3f &dir,     /// ray direction
    const Eigen::Vector3f &vmin,    /// lowest vertex
    const Eigen::Vector3f &vmax     /// highest vertex
);


DEFINE_TEST(ray_box_simple_intersection) {
    auto result = ray_intersects_box(
        Vec3(-2.0f, +0.0f, +0.0f),
        Vec3(+1.0f, +0.0f, +0.0f),
        BOX[0], BOX[1]
    );
    CHECK(true == result);
}


DEFINE_TEST(ray_box_no_intersection) {
    // Wrong ray direction
    auto result = ray_intersects_box(
        Vec3(-2.0f, +0.0f, +0.0f),
        Vec3(+0.0f, +1.0f, +0.0f),
        BOX[0], BOX[1]
    );
    CHECK(false == result);
}


//
// OBJ loading
//


DEFINE_TEST(load_wavefront) {
    std::istringstream ss(
        "v 0.0 0.0 0.0\n"
        "v 1.0 0.0 0.0\n"
        "v  0.0 1.0  0.0\n"
        "vn 1.0 1.0 0.0\n"
        "vn 0.0 1.0  0.0\n"
        "f 1/0/0 2/0/0  3/0/0\n"
        "f 2 3  1\n"
    );
    TriangleMesh mesh = TriangleMesh::from_wavefront_obj(ss);
    CHECK( mesh.vertices.size() == 3 );
    CHECK( mesh.vix_tuples.size() == 2 );
    CHECK( mesh.normals.size() == 2 );
}


//
// Binner3D
//


DEFINE_TEST(binner_3D_create) {
    using namespace Eigen;
    using Vector3u = Eigen::Matrix<std::size_t, 3, 1>;

    std::vector<Vector3f> vertices({
        Vector3f(0, 0, 0),
        Vector3f(2, 2, 2)
    });

    auto b0 = Binner3D::with_bin_size(1.f, vertices);
    CHECK( b0.bin_counts().isApprox(Vector3u(2, 2, 2)) );

    auto b1 = Binner3D::with_bin_count(8, vertices);
    CHECK( cppkit::mathx::is_close(b1.binsize, 1.f)  );
    CHECK( b1.bin_counts().isApprox(Vector3u(2, 2, 2)) );
}


DEFINE_TEST(binner_3D_binning) {
    using namespace Eigen;

    std::vector<Vector3f> vertices({
        Vector3f(0, 0, 0),
        Vector3f(2, 2, 2)
    });

    auto b0 = Binner3D::with_bin_size(1.f, vertices);
    CHECK( b0.bin(Vector3f(0.1f, 0.1f, 0.1f)).isApprox(Vector3i(0, 0, 0)) );
    CHECK( b0.bin(Vector3f(1.1f, 1.1f, 1.1f)).isApprox(Vector3i(1, 1, 1)) );
    CHECK( b0.bin(Vector3f(1.0f, 1.0f, 1.0f)).isApprox(Vector3i(1, 1, 1)) );
    CHECK( b0.bin(Vector3f(2.0f, 2.0f, 2.0f)).isApprox(Vector3i(2, 2, 2)) );

    auto d = b0.box_domain(Vector3f(0.1, 0.1, 0.1), Vector3f(1.5, 1.5, 1.5));
    CHECK( d.first.isApprox(Vector3i(0, 0, 0)) && d.second.isApprox(Vector3i(2, 2, 2)) );

    d = b0.box_domain(Vector3f(0.1, 0.1, 0.1), Vector3f(0.99, 0.99, 0.99));
    CHECK( d.first.isApprox(Vector3i(0, 0, 0)) && d.second.isApprox(Vector3i(1, 1, 1)) );

    d = b0.box_domain(Vector3f(0.1, 0.1, 0.1), Vector3f(0.9, 2.1, 3.0));
    CHECK( d.first.isApprox(Vector3i(0, 0, 0)) && d.second.isApprox(Vector3i(1, 3, 3)) );
}


//
// Ray casting
//


DEFINE_TEST(ray_cast_regression_test) {
    std::ifstream model_obj("models/__ray_cast_regression_test_0__.obj");
    IndexedTriangleMesh mesh = IndexedTriangleMesh::from_wavefront_obj(model_obj, 125);

    const auto &pos = Eigen::Vector3f(-1.00, 0.00, 0.25);
    const auto &orientation = Eigen::Vector3f(0, 0, 0);
    const auto &scan_dirs = LidarSensorScanGenerator::get_scan(pos, orientation);
    const auto MAX_RAY_DIST = 30.f;

    std::vector<Eigen::Vector3f> results;
    for (const auto &dir : scan_dirs) {
        const float dist = mesh.cast_ray(pos, dir);
        if (dist >= MAX_RAY_DIST) continue;
        results.push_back(pos + dist * dir);
    }

    std::ifstream result_obj("models/__ray_cast_regression_result_0__.obj");
    TriangleMesh expected = TriangleMesh::from_wavefront_obj(result_obj);

    CHECK( results.size() == expected.vertices.size()  );
}


int main() {
    puts("All tests passed");
    return 0;
}