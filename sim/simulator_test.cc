#include <cstdio>

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
        "v 0.0 1.0 0.0\n"
        "vn 1.0 1.0 0.0\n"
        "vn 0.0 1.0 0.0\n"
        "f 1/0/0 2/0/0 3/0/0\n"
        "f 2 3 1\n"
    );
    TriangleMesh mesh = TriangleMesh::from_wavefront_obj(ss);
    CHECK( mesh.vertices.size() == 3 );
    CHECK( mesh.indices.size() == 6 );
    CHECK( mesh.normals.size() == 2 );
}


int main() {
    puts("All tests passed");
    return 0;
}