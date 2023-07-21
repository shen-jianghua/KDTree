#pragma once

/*
 * file: KDTree.hpp
 * author: 沈江华
 *
 * This is an adaptation of the KD-tree implementation in rosetta code
 *  https://rosettacode.org/wiki/K-d_tree
 * It is a reimplementation of the C code using C++.
 * It also includes a few more queries than the original
 *
 */

#include <algorithm>
#include <functional>
#include <memory>
#include <vector>

using point_t = std::vector<double>;
using indexArr = std::vector<size_t>;
using pointIndex = typename std::pair<std::vector<double>, size_t>;

class KDNode 
{
public:
    using KDNodePtr = std::shared_ptr<KDNode>;
    size_t m_index;
    point_t m_x;
    KDNodePtr m_left;
    KDNodePtr m_right;

    // initializer
    KDNode();
    KDNode(const point_t&, const size_t&, const KDNodePtr&, const KDNodePtr&);
    KDNode(const pointIndex&, const KDNodePtr&, const KDNodePtr&);
    ~KDNode();

    // getter
    double coord(const size_t &);

    // conversions
    //点坐标不为空，返回True
    explicit operator bool();
    explicit operator point_t();
    explicit operator size_t();
    explicit operator pointIndex();
};

using KDNodePtr = std::shared_ptr<KDNode>;

KDNodePtr NewKDNodePtr();

// square euclidean distance
inline double dist2(const point_t&, const point_t&);
inline double dist2(const KDNodePtr&, const KDNodePtr&);

// euclidean distance
inline double dist(const point_t&, const point_t&);
inline double dist(const KDNodePtr&, const KDNodePtr&);

// Need for sorting
class comparer
{
public:
    size_t m_idx;//维度
    explicit comparer(size_t idx_);
    //按照第m_idx维度的坐标分量比较大小
    inline bool compare_idx(const pointIndex& a,
                            const pointIndex& b)
    {
        return (a.first.at(m_idx) < b.first.at(m_idx));
    }
};

using pointIndexArr = typename std::vector< pointIndex >;

//按照第idx维度的坐标分量排序，从小到大
inline void sort_on_idx(
    const pointIndexArr::iterator&,
    const pointIndexArr::iterator&,
    size_t idx);

using pointVec = std::vector<point_t>;

class KDTree
{
private:
    KDNodePtr m_root;
    KDNodePtr m_leaf;

public:
    KDTree() = default;
    explicit KDTree(pointVec point_array);//explicit关键字的作用就是防止类构造函数的隐式自动转换

public:
    //查找最近点
    point_t nearest_point(const point_t& pt);
    //查找最近点的索引
    size_t nearest_index(const point_t& pt);
    //查找最近点和索引
    pointIndex nearest_pointIndex(const point_t& pt);

    //查找邻域点和索引
    pointIndexArr neighborhood(const point_t& pt, const double& rad);
    //查找邻域点
    pointVec neighborhood_points(const point_t& pt, const double& rad);
    //查找邻域点的索引
    indexArr neighborhood_indices(const point_t& pt, const double& rad);

private:
    KDNodePtr make_tree(
        const pointIndexArr::iterator& begin,
        const pointIndexArr::iterator& end,
        const size_t& length,
        const size_t& level);
    /// <summary>
    /// 查找最近点
    /// </summary>
    /// <param name="branch">       KDtree的分支(从root开始)</param>
    /// <param name="pt">           查询点</param>
    /// <param name="level">        从第level层开始(即在第level维度分割左右子树的)</param>
    /// <param name="best">         最近的点(从根节点开始)</param>
    /// <param name="best_dist">    最近的距离</param>
    /// <returns></returns>
    KDNodePtr nearest_(
        const KDNodePtr& branch,
        const point_t& pt,
        const size_t& level,
        const KDNodePtr& best,
        const double& best_dist);

    // default caller
    KDNodePtr nearest_(const point_t& pt);

    pointIndexArr neighborhood_(
        const KDNodePtr& branch,
        const point_t& pt,
        const double& rad,
        const size_t& level);
};
