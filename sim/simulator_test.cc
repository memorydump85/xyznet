#include <cstdio>
#include <fstream>

#include <cppkit/mathx.hh>
#include <cppkit/test_framework.hh>

#include "simulator.hh"



using Vec3 = Eigen::Vector3f;

const Vec3 TRIANGLE[3] = {
    Vec3(+0.0f, +0.0f, +0.0f),
    Vec3(+1.0f, +0.0f, +0.0f),
    Vec3(+0.0f, +1.0f, +0.0f),
};


DEFINE_TEST(simple_intersection) {
    auto dist = intersect_ray_with_triangle(
        Vec3(+0.5f, +0.5f, +1.0f),
        Vec3(+0.0f, +0.0f, -1.0f),
        TRIANGLE[0], TRIANGLE[1], TRIANGLE[2]
    );
    CHECK(cppkit::mathx::is_close(dist, 1.0f));
}


#include <iostream>
DEFINE_TEST(no_intersection) {
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
    CHECK( mesh.indices.size() == 6 );
    CHECK( mesh.normals.size() == 2 );
}


DEFINE_TEST(ray_cast_regression_test) {
    std::ifstream model_obj("models/__ray_cast_regression_test_0__.obj");
    TriangleMesh mesh = TriangleMesh::from_wavefront_obj(model_obj);

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