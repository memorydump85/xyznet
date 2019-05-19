#include <limits>
#include <unordered_set>
#include <unordered_map>

#include <cppkit/validation_macros.h>
#include <cppkit/string.hh>
#include <cppkit/mathx.hh>
#include <cppkit/rigidbody.hh>

#include "simulator.hh"



//
// Function implementations
//


/// Intersect `ray` with `triangle`
///
/// Returns the distance of the point of intersection from `origin`;
/// Returns `std::numeric_limits<float>::infinity()` if there is no
/// intersection.
///
/// NOTE: Implementation uses the Möller–Trumbore intersection
/// algorithm.
///
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


/// Does ray intersect (axis-aligned) box?
///
bool ray_intersects_box(
    const Eigen::Vector3f &origin,  /// ray origin
    const Eigen::Vector3f &dir,     /// ray direction
    const Eigen::Vector3f &vmin,    /// lowest vertex
    const Eigen::Vector3f &vmax     /// highest vertex
)
{
#ifdef UNROLL_RAY_BOX_INTERSECTION
    const float tx0 = (vmin.x() - origin.x()) / dir.x();
    const float tx1 = (vmax.x() - origin.x()) / dir.x();
    const float ty0 = (vmin.y() - origin.y()) / dir.y();
    const float ty1 = (vmax.y() - origin.y()) / dir.y();
    const float tz0 = (vmin.z() - origin.z()) / dir.z();
    const float tz1 = (vmax.z() - origin.z()) / dir.z();

    using std::min;
    using std::max;
    const float t_max = max({ min(tx0, tx1), min(ty0, ty1), min(tz0, tz1) });
    const float t_min = min({ max(tx0, tx1), max(ty0, ty1), max(tz0, tz1) });
    return (t_min > 0) && (t_max < t_min);
#else
    //
    // This implementation is described in:
    //
    // A Ray-Box Intersection Algorithm and Efficient Dynamic Voxel Rendering
    // by Majercik, Crassin, Shirley, and McGuire
    //
    const auto &inv_dir = dir.array().inverse();
    const auto &t0 = (vmin - origin).array() * inv_dir;
    const auto &t1 = (vmax - origin).array() * inv_dir;
    const auto &tmin = t0.cwiseMin(t1);
    const auto &tmax = t0.cwiseMax(t1);

    return tmin.maxCoeff() <= tmax.minCoeff();
#endif
}


//
// TriangleMesh method implementations
//


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


//
// LidarSensorScanGenerator method implementations
//


const float LidarSensorScanGenerator::START_ANGLE = cppkit::mathx::to_radians(+185.f);
const float LidarSensorScanGenerator::ANGLE_RANGE = cppkit::mathx::to_radians(350.f);


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


//
// XYZGrid method implementations
//


XYZGrid::XYZGrid(
        const Eigen::Vector3f &lbound_,
        const Eigen::Vector3f &ubound_,
        const float &cellsize_)
    : lbound(lbound_)
    , ubound(ubound_)
    , cellsize(cellsize_)
{
    CHECK( cellsize > 0 );

    const auto &cell_counts = ((ubound - lbound) / cellsize).array().floor().cast<int>();
    CHECK_MSG( cell_counts.minCoeff() >= 0,
        cppkit::strfmt("counts = {%d, %d, %d}", cell_counts(0), cell_counts(1), cell_counts(2)) );
    CHECK_MSG( cell_counts.minCoeff() < 0xffff,
        cppkit::strfmt("`cellsize == %f` is too small", cellsize) );
}


template<class C>
static inline auto bounds_(const C& vertices) {
    Eigen::Vector3f lbound = *(vertices.cbegin());
    Eigen::Vector3f ubound = *(vertices.cbegin());
    for (const auto &v : vertices) {
        lbound = lbound.cwiseMin(v);
        ubound = ubound.cwiseMax(v);
    }
    return std::make_pair(lbound, ubound);
}


template<class C>
XYZGrid XYZGrid::with_cell_size(
    const float &cellsize,
    const C& vertices )
{
    const auto &[lbound, ubound] = bounds_(vertices);
    return XYZGrid(lbound, ubound, cellsize);
}


template<class C>
XYZGrid XYZGrid::with_cell_count(
    const std::size_t &cell_count,
    const C& vertices)
{
    const auto &[lbound, ubound] = bounds_(vertices);
    const auto& cellsize = pow((ubound - lbound).prod() / cell_count, 1./3.);
    return XYZGrid(lbound, ubound, cellsize);
}


//
// IndexedTriangleMesh method implementations
//


IndexedTriangleMesh::IndexedTriangleMesh(
    const std::vector<Eigen::Vector3f> &vx,
    const std::vector<Eigen::Vector3i> &ix,
    const std::vector<Eigen::Vector3f> &nm,
    const std::size_t &cell_count )
    : vertices(vx)
    , normals(nm)
    , ix_tuples(ix)
    , grid_(XYZGrid::with_cell_count(cell_count, vx))
{
    for (const auto &ix : ix_tuples) {
        const auto &v0 = vertices[ix(0)];
        const auto &v1 = vertices[ix(1)];
        const auto &v2 = vertices[ix(2)];

        // Add triangle to covered cells
        const auto &vmin = v0.cwiseMin(v1).cwiseMin(v2);
        const auto &vmax = v0.cwiseMax(v1).cwiseMax(v2);

        const auto &[r, s] = grid_.covered_cell_range(vmin, vmax);
        CHECK( r.minCoeff() >= 0 );
        CHECK( s.minCoeff() >= 0 );

        for (int i = r.x(); i < s.x(); ++i)
            for (int j = r.y(); j < s.y(); ++j)
                for (int k = r.z(); k < s.z(); ++k) {
                    cells_[grid_.flat_ix(i, j, k)].push_back(ix);
                }
    }
}


IndexedTriangleMesh IndexedTriangleMesh::from_wavefront_obj(
    std::istream& is,
    const std::size_t &cell_count)
{
    TriangleMesh mesh = TriangleMesh::from_wavefront_obj(is);
    return IndexedTriangleMesh(mesh.vertices, mesh.ix_tuples, mesh.normals, cell_count);
}


float IndexedTriangleMesh::cast_ray(
    const Eigen::Vector3f &origin,
const Eigen::Vector3f &dir ) const
{
    const auto &INF = std::numeric_limits<float>::infinity();

    if (false == ray_intersects_box(origin, dir, grid_.lbound, grid_.ubound)) {
        return INF;
    }

    float min_dist = INF;
    for (const auto &c : cells_) {
        const auto &[cmin, cmax] = grid_.cell_bounds(grid_.lift_ix(c.first));
        if (false == ray_intersects_box(origin, dir, cmin, cmax)) {
            continue;
        }

        for (const auto &ix : c.second) {
            const auto& v0 = vertices[ix(0)];
            const auto& v1 = vertices[ix(1)];
            const auto& v2 = vertices[ix(2)];
            const float &dist = intersect_ray_with_triangle(origin, dir, v0, v1, v2);
            min_dist = std::min(min_dist, dist);
        }
    }

    return min_dist;
}
