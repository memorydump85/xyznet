#pragma once

#include <limits>
#include <unordered_set>

#include <cppkit/string.hh>
#include <Eigen/Dense>



/// Intersect `ray` with `triangle`
///
/// Returns the distance of the point of intersection from `origin`;
/// Returns `std::numeric_limits<float>::infinity()` if there is no
/// intersection.
///
/// NOTE: Implementation uses the Möller–Trumbore intersection
/// algorithm.
///
static inline
float intersect_ray_with_triangle(
    const Eigen::Vector3f &origin,  /// ray origin
    const Eigen::Vector3f &dir,     /// ray direction
    const Eigen::Vector3f &v0,      /// triangle vertex1
    const Eigen::Vector3f &v1,      /// triangle vertex2
    const Eigen::Vector3f &v2       /// triangle vertex3
);


//-------------------------------------
class TriangleMesh {
//-------------------------------------
public:
    const std::vector<Eigen::Vector3f> vertices;
    const std::vector<Eigen::Vector3f> normals;
    const std::vector<std::size_t> indices;

public:
    TriangleMesh(
        const std::vector<Eigen::Vector3f> &vx,
        const std::vector<std::size_t> &ix,
        const std::vector<Eigen::Vector3f> &nm )
        : vertices(vx)
        , normals(nm)
        , indices(ix)
        { CHECK( (indices.size() % 3) == 0 ); }

    /// Load mesh from .OBJ file
    static TriangleMesh from_wavefront_obj(std::istream& is);

    /// Perform ray casting
    /// return distance of intersection from ray origin
    float cast_ray(const Eigen::Vector3f &origin, const Eigen::Vector3f &dir) const;

    /// get point of intersection of ray with mesh
    Eigen::Vector3f intersect_ray(const Eigen::Vector3f &origin, const Eigen::Vector3f &dir ) const
        { return origin + dir * cast_ray(origin, dir); }
};



//
// Function implementations
//

#include <iostream>
float intersect_ray_with_triangle(
    const Eigen::Vector3f &origin,  /// ray origin
    const Eigen::Vector3f &dir,     /// ray direction
    const Eigen::Vector3f &v0,      /// triangle vertex1
    const Eigen::Vector3f &v1,      /// triangle vertex2
    const Eigen::Vector3f &v2       /// triangle vertex3
)
{
    const float EPSILON = 1e-8f;
    const auto& INF = std::numeric_limits<float>::infinity();

    const auto &v0v1 = v1 - v0;
    const auto &v0v2 = v2 - v0;
    const auto &pvec = dir.cross(v0v2);
    const auto &det = v0v1.dot(pvec);

    if (fabs(det) < EPSILON) /* Ray is parallel to triangle and misses it */ {
        return INF;
    }

    const float &invDet = 1 / det;

    const auto &tvec = origin - v0;
    const auto &u = tvec.dot(pvec) * invDet;
    if (u < 0 || u > 1) return INF;

    const auto &qvec = tvec.cross(v0v1);
    const auto &v = dir.dot(qvec) * invDet;
    if (v < 0 || u + v > 1) return INF;

    const auto &t = v0v2.dot(qvec) * invDet;
    return t > EPSILON ? t : INF;
}


//
// NOTE: This implementation may not be robust when presented with a
// malformed OBJ source.
//
TriangleMesh TriangleMesh::from_wavefront_obj(std::istream& is) {
    std::size_t linenum = 0;
    std::string line;

    enum {
        ERR_INVALID_LINE    = 1,
        ERR_VP_UNSUPPORTED  = 2,
        ERR_TRIFACES_ONLY   = 3,
        ERR_UNEXPECTED      = 4,
    };

    using cppkit::strfmt;
    const auto& fail = [&] (const auto &ecode) {
        const auto &context = strfmt("line %d -- ", linenum);
        switch(ecode) {
        case ERR_INVALID_LINE:
            CHECK_MSG(false, context + "cannot parse line: '" + line + "'");
            break;
        case ERR_VP_UNSUPPORTED:
            CHECK_MSG(false, context + "parameter-space vertices are not supported");
            break;
        case ERR_TRIFACES_ONLY:
            CHECK_MSG(false, context + "only triangular faces are supported");
            break;
        case ERR_UNEXPECTED:
            CHECK_MSG(false, context + "unexpected line");
            break;
        default:
            CHECK(false);
        }
    };

    const std::unordered_set<std::string> ignored {
        "o", "g", "s", "vt", "vp", "l", "mtllib", "usemtl"
    };

    std::vector<Eigen::Vector3f> vertices;
    std::vector<std::size_t> indices;
    std::vector<Eigen::Vector3f> normals;

    const auto part1 = [] (const std::string &s, const char sep = '/') {
        return s.substr(0, std::min(s.size(), s.find(sep)));
    };

    while(std::getline(is, line)) {
        ++linenum;
        if(line.empty() || line[0] == '#') continue;

        cppkit::StringTokenGenerator tokens(line, ' ');
        if (!tokens.has_next()) fail( ERR_INVALID_LINE );

        const auto &tok0 = tokens.next();
        if (tok0 == "v") {
            float x = std::stof(tokens.next());
            float y = std::stof(tokens.next());
            float z = std::stof(tokens.next());
            vertices.emplace_back(x, y, z);
        }
        else if (tok0 == "vn") {
            float i = std::stof(tokens.next());
            float j = std::stof(tokens.next());
            float k = std::stof(tokens.next());
            normals.emplace_back(i, j, k);
        }
        else if (tok0 == "f") {
            // NOTE: OBJ file-format specifies 1-based vertex indices
            std::size_t i = std::stoul(part1(tokens.next())) - 1;
            CHECK( i < vertices.size() );
            indices.push_back(i);

            std::size_t j = std::stoul(part1(tokens.next())) - 1;
            CHECK( j < vertices.size() );
            indices.push_back(j);

            std::size_t k = std::stoul(part1(tokens.next())) - 1;
            CHECK( k < vertices.size() );
            indices.push_back(k);

            if (tokens.has_next()) fail( ERR_TRIFACES_ONLY );
        }
        else if (tok0 == "vp") {
            fail ( ERR_VP_UNSUPPORTED );
        }
        else if (0 != ignored.count(tok0)) {
            continue;
        }
        else {
            fail( ERR_UNEXPECTED );
        }
    }

    return std::move(TriangleMesh(vertices, indices, normals));
}


float TriangleMesh::cast_ray(
    const Eigen::Vector3f &origin,
    const Eigen::Vector3f &dir
) const
{
    float min_dist = std::numeric_limits<float>::infinity();
    for (std::size_t i=0; i < indices.size(); i+=3) {
        const auto& v0 = vertices[indices[i]];
        const auto& v1 = vertices[indices[i+1]];
        const auto& v2 = vertices[indices[i+2]];
        const float &dist = intersect_ray_with_triangle(origin, dir, v0, v1, v2);
        min_dist = std::min(min_dist, dist);
    }
    return min_dist;
}
