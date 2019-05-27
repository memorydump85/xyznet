#pragma once

#include <Eigen/Dense>



//-------------------------------------
class TriangleMesh {
//-------------------------------------
public:
    const std::vector<Eigen::Vector3f> vertices;
    const std::vector<Eigen::Vector3f> normals;
    const std::vector<Eigen::Vector3i> vix_tuples;

public:
    TriangleMesh(
        const std::vector<Eigen::Vector3f> &vx,
        const std::vector<Eigen::Vector3i> &ix,
        const std::vector<Eigen::Vector3f> &nm )
        : vertices(vx)
        , normals(nm)
        , vix_tuples(ix)
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


///
/// Quantize 3 dimensional values into bins of a specified size.
///
//------------------------------------
class Binner3D {
//-------------------------------------
    using SizeVector3 = Eigen::Matrix<std::size_t, 3, 1>;

public:
    const Eigen::Vector3f lbound;   /// [ min(x), min(y), min(z) ]
    const Eigen::Vector3f ubound;   /// [ max(x), max(y), max(z) ]
    const float binsize;            /// resolution. same for all dimensions

    Binner3D(
        const Eigen::Vector3f &lbound,
        const Eigen::Vector3f &ubound,
        const float &binsize );

    static Binner3D with_bin_size(
        const float &binsize,
        const std::vector<Eigen::Vector3f>& vertices );

    /// Create a Binner3D object which bins values into `bin_count`
    /// total bins. Bins are constructed such that each dimension has
    /// `pow(bin_count, 1.0/3)` bins.
    static Binner3D with_bin_count(
        const std::size_t &bin_count,
        const std::vector<Eigen::Vector3f>& vertices );


    /// Get bin index for point
    Eigen::Vector3i bin(const Eigen::Vector3f &p) const
        { return ((p - lbound) / binsize).array().floor().cast<int>(); }

    /// Get the range of values that are discretized into the bin
    /// represented by (i, j, k)
    std::pair<Eigen::Vector3f, Eigen::Vector3f>
    bin_range(int i, int j, int k) const {
        return std::make_pair(
            lbound + Eigen::Vector3f(i,   j,   k  ) * binsize,
            lbound + Eigen::Vector3f(i+1, j+1, k+1) * binsize
        );
    }

    /// Get the range of values that are discretized into the bin
    /// represented by (i, j, k)
    std::pair<Eigen::Vector3f, Eigen::Vector3f>
    bin_range(const Eigen::Vector3i &i) const
        { return bin_range(i(0), i(1), i(2)); }

    /// Get bins covered by the box with inclusive lower bound `vmin` and
    /// exclusive upper bound `vmax`.
    std::pair<Eigen::Vector3i, Eigen::Vector3i>
    box_domain(
        const Eigen::Vector3f &vmin,
        const Eigen::Vector3f &vmax ) const
    {
        const auto &lo = ((vmin - lbound) / binsize).array().floor().cast<int>();
        const auto &hi = ((vmax - lbound) / binsize).array().ceil().cast<int>();
        // `hi` should be at least 1 more than `lo` in all dimensions
        return std::make_pair(lo, lo + (hi - lo).cwiseMax(1));
    }

    /// Get dimensional bin counts
    SizeVector3 bin_counts() const {
        return ((ubound - lbound) / binsize).array().ceil().cast<std::size_t>().cwiseMax(1);
    }
};


///
/// Converter that translates between 3D and linear indices.
///
/// A typical implementation would use a 3 dimensional size
/// specifications to translate from 3D to linear indices. In this
/// implementation, however, we assume that the 3D index values can be
/// represented using at most 20 bits and hence are in the closed
/// interval [0, 2^21 - 1]
///
//-------------------------------------
class IndexMapper3D {
//-------------------------------------
public:
    static const std::size_t MAX_INDEX = 0xfffff;

    /// Convert 3D index to linear index
    std::size_t flatten(std::size_t i, std::size_t j, std::size_t k) const {
        CHECK( (i <= MAX_INDEX) && (j <= MAX_INDEX) && (k <= MAX_INDEX) );
        { return (i << 40) | (j << 20) | k; }
    }

    /// Convert 3D index to linear index
    uint64_t flatten(const Eigen::Vector3i &v) const {
        CHECK( v.minCoeff() >= 0 );
        return flatten(v(0), v(1), v(2));
    }

    /// Extract 3D index from linear index
    Eigen::Vector3i unravel(const uint64_t &i) const {
        const std::size_t MASK_20LSB = 0xfffff;
        return Eigen::Vector3i((i >> 40) & MASK_20LSB, (i >> 20) & MASK_20LSB, i & MASK_20LSB);
    }
};


//-------------------------------------
class IndexedTriangleMesh {
//-------------------------------------
public:
    const std::vector<Eigen::Vector3f> vertices;
    const std::vector<Eigen::Vector3f> normals;
    const std::vector<Eigen::Vector3i> vix_tuples;  /// Triangle vertex indices

private:
    const Binner3D binner_;
    const IndexMapper3D ix_mapper_;
    using bin_contents_t = std::vector<Eigen::Vector3i>;
    std::unordered_map<uint64_t, bin_contents_t> bins_;

public:
    IndexedTriangleMesh(
        const std::vector<Eigen::Vector3f> &vx,
        const std::vector<Eigen::Vector3i> &ix,
        const std::vector<Eigen::Vector3f> &nm,
        const std::size_t &bin_count );

    /// Load mesh from .OBJ file
    static
    IndexedTriangleMesh from_wavefront_obj(
        std::istream& is,
        const std::size_t &bin_count );

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
