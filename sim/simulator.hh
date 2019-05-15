#pragma once

#include <limits>
#include <unordered_set>

#include <cppkit/string.hh>
#include <cppkit/mathx.hh>
#include <cppkit/rigidbody.hh>
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
    const std::vector<Eigen::Vector3i> ix_tuples;

public:
    TriangleMesh(
        const std::vector<Eigen::Vector3f> &vx,
        const std::vector<Eigen::Vector3i> &ix,
        const std::vector<Eigen::Vector3f> &nm )
        : vertices(vx)
        , normals(nm)
        , ix_tuples(ix)
        { }

    /// Load mesh from .OBJ file
    static TriangleMesh from_wavefront_obj(std::istream& is);

    /// Perform ray casting
    /// return distance of intersection from ray origin
    float cast_ray(const Eigen::Vector3f &origin, const Eigen::Vector3f &dir) const;

    /// get point of intersection of ray with mesh
    Eigen::Vector3f intersect_ray(const Eigen::Vector3f &origin, const Eigen::Vector3f &dir ) const
        { return origin + dir * cast_ray(origin, dir); }
};


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
        , rot_(cppkit::mathx::rotation_matrix(rph(0), rph(1), rph(2)).block(0, 0, 3, 3))
    {}

    bool has_next() const
        { return ray_ix_ < (NUM_SCANS * NUM_RAYS_PER_SCAN); }

    /// Peek the next ray direction
    Eigen::Vector3f peek() const;

    /// Return the next ray direction
    Eigen::Vector3f next();

    /// Get a vector containing rays for an entire scan
    static
    std::vector<Eigen::Vector3f> get_scan(
        const Eigen::Vector3f &pos,
        const Eigen::Vector3f &orientation );
};

const float LidarSensorScanGenerator::START_ANGLE = cppkit::mathx::to_radians(+185.f);
const float LidarSensorScanGenerator::ANGLE_RANGE = cppkit::mathx::to_radians(350.f);



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
    std::vector<Eigen::Vector3i> ix_tuples;
    std::vector<Eigen::Vector3f> normals;

    const auto part1 = [] (const std::string &s, const char sep = '/') {
        return s.substr(0, std::min(s.size(), s.find(sep)));
    };

    while(std::getline(is, line)) {
        ++linenum;
        if(line.empty() || line[0] == '#') continue;

        cppkit::StringTokenGenerator tokens(line, ' ', /*ignore_consecutive*/ true);
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
            // NOTE: OBJ file-format specifies 1-based vertex ix_tuples
            std::size_t i = std::stoul(part1(tokens.next())) - 1;
            std::size_t j = std::stoul(part1(tokens.next())) - 1;
            std::size_t k = std::stoul(part1(tokens.next())) - 1;

            CHECK( i < vertices.size() );
            CHECK( j < vertices.size() );
            CHECK( k < vertices.size() );
            ix_tuples.emplace_back(i, j, k);

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

    return std::move(TriangleMesh(vertices, ix_tuples, normals));
}


float TriangleMesh::cast_ray(
    const Eigen::Vector3f &origin,
    const Eigen::Vector3f &dir
) const
{
    float min_dist = std::numeric_limits<float>::infinity();
    for (const auto &ix : ix_tuples) {
        const auto& v0 = vertices[ix(0)];
        const auto& v1 = vertices[ix(1)];
        const auto& v2 = vertices[ix(2)];
        const float &dist = intersect_ray_with_triangle(origin, dir, v0, v1, v2);
        min_dist = std::min(min_dist, dist);
    }
    return min_dist;
}


Eigen::Vector3f LidarSensorScanGenerator::peek() const {
    const auto &DEL_THETA = ANGLE_RANGE / NUM_RAYS_PER_SCAN;
    const auto &theta = START_ANGLE + (ray_ix_ % NUM_RAYS_PER_SCAN) * DEL_THETA;

    const auto &scan_ix = ray_ix_ / NUM_RAYS_PER_SCAN;
    const auto &phi = cppkit::mathx::to_radians(90 + 3*(-16.f + scan_ix));

    const auto &xyz = Eigen::Vector3f(cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi));
    return rot_ * xyz;
}


Eigen::Vector3f LidarSensorScanGenerator::next() {
    const auto& dir = peek();
    ++ray_ix_;
    return dir;
}


std::vector<Eigen::Vector3f> LidarSensorScanGenerator::get_scan(
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
