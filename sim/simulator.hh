#pragma once

#include <Eigen/Dense>



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
    float cast_ray(
        const Eigen::Vector3f &origin,
        const Eigen::Vector3f &dir ) const;

    /// get point of intersection of ray with mesh
    Eigen::Vector3f intersect_ray(
        const Eigen::Vector3f &origin,
        const Eigen::Vector3f &dir ) const
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


//-------------------------------------
class XYZGrid {
//-------------------------------------
public:
    const Eigen::Vector3f lbound;  /// [ min(x), min(y), min(z) ]
    const Eigen::Vector3f ubound;  /// [ max(x), max(y), max(z) ]
    const float cellsize;

    XYZGrid(
        const Eigen::Vector3f &lbound,
        const Eigen::Vector3f &ubound,
        const float &cellsize );

    template<class C>
    static XYZGrid with_cell_size(
        const float &cellsize,
        const C& vertices );

    template<class C>
    static XYZGrid with_cell_count(
        const std::size_t &cell_count,
        const C& vertices );

    /// XYZGrid index to linear index
    uint64_t flat_ix(uint64_t i, uint64_t j, uint64_t k) const
        { return (i << 32) | (j << 16) | k; }

    /// XYZGrid index to linear index
    uint64_t flat_ix(const Eigen::Vector3i &v) const
        { return flat_ix(v(0), v(1), v(2)); }

    /// Extract grid index from linear index
    Eigen::Vector3i lift_ix(const uint64_t &i) const
        { return Eigen::Vector3i((i >> 32) & 0xffff, (i >> 16) & 0xffff, i & 0xffff); }

    /// Get bounds of the grid cell at (i, j, k)
    std::pair<Eigen::Vector3f, Eigen::Vector3f>
    cell_bounds(int i, int j, int k) const {
        return std::make_pair(
            lbound + Eigen::Vector3f(i,   j,   k  ) * cellsize,
            lbound + Eigen::Vector3f(i+1, j+1, k+1) * cellsize
        );
    }

    /// Get bounds of the grid cell at (i, j, k)
    std::pair<Eigen::Vector3f, Eigen::Vector3f>
    cell_bounds(const Eigen::Vector3i &i) const
        { return cell_bounds(i(0), i(1), i(2)); }

    /// Get cells covered by the box with lower bound `vmin` and
    /// upper bound `vmax`.
    std::pair<Eigen::Vector3i, Eigen::Vector3i>
    covered_cell_range(
        const Eigen::Vector3f &vmin,
        const Eigen::Vector3f &vmax ) const
    {
        return std::make_pair(
            ((vmin - lbound) / cellsize).array().floor().cast<int>(),
            ((vmax - lbound) / cellsize).array().ceil().cast<int>() + 1
        );
    }
};


//-------------------------------------
class IndexedTriangleMesh {
//-------------------------------------
public:
    const std::vector<Eigen::Vector3f> vertices;
    const std::vector<Eigen::Vector3f> normals;
    const std::vector<Eigen::Vector3i> ix_tuples;

private:
    const XYZGrid grid_;
    using cell_contents_t_ = std::vector<Eigen::Vector3i>;
    std::unordered_map<uint64_t, cell_contents_t_> cells_;

public:
    IndexedTriangleMesh(
        const std::vector<Eigen::Vector3f> &vx,
        const std::vector<Eigen::Vector3i> &ix,
        const std::vector<Eigen::Vector3f> &nm,
        const std::size_t &cell_count );

    /// Load mesh from .OBJ file
    static
    IndexedTriangleMesh from_wavefront_obj(
        std::istream& is,
        const std::size_t &cell_count );

    /// Perform ray casting
    /// return distance of intersection from ray origin
    float cast_ray(
        const Eigen::Vector3f &origin,
        const Eigen::Vector3f &dir ) const;

    /// get point of intersection of ray with mesh
    Eigen::Vector3f intersect_ray(
        const Eigen::Vector3f &origin,
        const Eigen::Vector3f &dir ) const
    { return origin + dir * cast_ray(origin, dir); }
};
