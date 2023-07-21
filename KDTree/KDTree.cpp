/*
 * file: KDTree.hpp
 * author: 沈江华
 *
 * This is an adaptation of the KD-tree implementation in rosetta code
 * https://rosettacode.org/wiki/K-d_tree
 *
 * It is a reimplementation of the C code using C++.  It also includes a few
 * more queries than the original, namely finding all points at a distance
 * smaller than some given distance to a point.
 *
 */

#include <algorithm>
#include <cmath>
#include <functional>
#include <iterator>
#include <limits>
#include <memory>
#include <vector>

#include "KDTree.hpp"

KDNode::KDNode() = default;

KDNode::KDNode(const point_t &pt, const size_t &idx_, const KDNodePtr &left_,
               const KDNodePtr &right_)
{
    m_x = pt;
    m_index = idx_;
    m_left = left_;
    m_right = right_;
}

KDNode::KDNode(const pointIndex &pi, const KDNodePtr &left_,
               const KDNodePtr &right_) 
{
    m_x = pi.first;
    m_index = pi.second;
    m_left = left_;
    m_right = right_;
}

KDNode::~KDNode() = default;

double KDNode::coord(const size_t &idx)
{ 
    return m_x.at(idx);
}

KDNode::operator bool()
{ 
    return (!m_x.empty());
}

KDNode::operator point_t() 
{ 
    return m_x;
}

KDNode::operator size_t() 
{ 
    return m_index;
}

KDNode::operator pointIndex()
{
    return pointIndex(m_x, m_index);
}

KDNodePtr NewKDNodePtr() 
{
    KDNodePtr mynode = std::make_shared< KDNode >();
    return mynode;
}

inline double dist2(const point_t &a, const point_t &b)
{
    double distc = 0;
    for (size_t i = 0; i < a.size(); i++)
    {
        double di = a.at(i) - b.at(i);
        distc += di * di;
    }
    return distc;
}

inline double dist2(const KDNodePtr &a, const KDNodePtr &b)
{
    return dist2(a->m_x, b->m_x);
}

inline double dist(const point_t &a, const point_t &b)
{
    return std::sqrt(dist2(a, b));
}

inline double dist(const KDNodePtr &a, const KDNodePtr &b)
{
    return std::sqrt(dist2(a, b));
}


comparer::comparer(size_t idx_)
    : m_idx{idx_}
{
}


inline void sort_on_idx(const pointIndexArr::iterator &begin,
                        const pointIndexArr::iterator &end,
                        size_t idx)
{
    comparer comp(idx);
    comp.m_idx = idx;

    using std::placeholders::_1;//占位符
    using std::placeholders::_2;
    //std::nth_element部分排序算法，它重排 [first, last) 中元素，第nth个元素之前的元素都小于第nth个元素，而第n个元素之后的元素大于第nth个元素
    //std::bind可调用对象与参数一起绑定为另一个std::function供调用
    std::nth_element(begin, begin + std::distance(begin, end) / 2,
                     end, std::bind(&comparer::compare_idx, comp, _1, _2));//bind to a pointer to member function
}

using pointVec = std::vector< point_t >;

KDTree::KDTree(pointVec point_array)
{
    m_leaf = std::make_shared<KDNode>();
    // iterators
    pointIndexArr arr;
    for (size_t i = 0; i < point_array.size(); i++)
    {
        arr.push_back(pointIndex(point_array.at(i), i));
    }

    auto begin = arr.begin();
    auto end = arr.end();

    size_t length = arr.size();
    size_t level = 0;  // starting

    m_root = KDTree::make_tree(begin, end, length, level);
}

KDNodePtr KDTree::make_tree(const pointIndexArr::iterator& begin,
                            const pointIndexArr::iterator& end,
                            const size_t& length,
                            const size_t& level)
{
    if (begin == end)
    {
        return NewKDNodePtr();  // empty tree
    }

    size_t dim = begin->first.size();

    if (length > 1) 
    {
        sort_on_idx(begin, end, level);
    }

    auto middle = begin + (length / 2);

    auto l_begin = begin;
    auto l_end = middle;
    auto r_begin = middle + 1;
    auto r_end = end;

    size_t l_len = length / 2;
    size_t r_len = length - l_len - 1;

    KDNodePtr left;
    if (l_len > 0 && dim > 0)
    {
        left = make_tree(l_begin, l_end, l_len, (level + 1) % dim);
    }
    else
    {
        left = m_leaf;
    }
    KDNodePtr right;
    if (r_len > 0 && dim > 0)
    {
        right = make_tree(r_begin, r_end, r_len, (level + 1) % dim);
    }
    else
    {
        right = m_leaf;
    }

    // KDNode result = KDNode();
    return std::make_shared< KDNode >(*middle, left, right);
}

KDNodePtr KDTree::nearest_(
    const KDNodePtr& branch,
    const point_t& pt,
    const size_t& level,
    const KDNodePtr& best,
    const double& best_dist)
{
    double d, dx, dx2;

    if (!bool(*branch))//判断节点是否为空
    {
        return NewKDNodePtr();  // basically, null
    }

    point_t branch_pt(*branch);//
    size_t dim = branch_pt.size();

    d = dist2(branch_pt, pt);
    dx = branch_pt.at(level) - pt.at(level);
    dx2 = dx * dx;

    KDNodePtr best_l = best;
    double best_dist_l = best_dist;

    if (d < best_dist) {
        best_dist_l = d;
        best_l = branch;
    }

    size_t next_lv = (level + 1) % dim;
    KDNodePtr section;
    KDNodePtr other;

    // select which branch makes sense to check
    if (dx > 0) {
        section = branch->m_left;
        other = branch->m_right;
    } else {
        section = branch->m_right;
        other = branch->m_left;
    }

    // keep nearest neighbor from further down the tree
    // 找到叶子节点
    KDNodePtr further = nearest_(section, pt, next_lv, best_l, best_dist_l);
    if (!further->m_x.empty()) {
        double dl = dist2(further->m_x, pt);
        if (dl < best_dist_l) {
            best_dist_l = dl;
            best_l = further;
        }
    }
    // only check the other branch if it makes sense to do so
    // 进行'回溯'操作：算法沿搜索路径反向查找是否有距离查询点更近的数据点
    if (dx2 < best_dist_l) {
        further = nearest_(other, pt, next_lv, best_l, best_dist_l);
        if (!further->m_x.empty()) {
            double dl = dist2(further->m_x, pt);
            if (dl < best_dist_l) {
                best_dist_l = dl;
                best_l = further;
            }
        }
    }

    return best_l;
}

// default caller
KDNodePtr KDTree::nearest_(const point_t &pt)
{
    size_t level = 0;
    // KDNodePtr best = branch;
    double branch_dist = dist2(point_t(*m_root), pt);
    return nearest_(m_root,          // beginning of tree
                    pt,            // point we are querying
                    level,         // start from level 0
                    m_root,          // best is the root
                    branch_dist);  // best_dist = branch_dist
}

pointIndexArr KDTree::neighborhood_(
    const KDNodePtr& branch,
    const point_t& pt,
    const double& rad,
    const size_t& level)
{
    double d, dx, dx2;

    if (!bool(*branch))
    {
        // branch has no point, means it is a leaf,
        // no points to add
        return pointIndexArr();
    }

    size_t dim = pt.size();

    double r2 = rad * rad;

    d = dist2(point_t(*branch), pt);
    dx = point_t(*branch).at(level) - pt.at(level);
    dx2 = dx * dx;

    pointIndexArr nbh, nbh_s, nbh_o;
    if (d <= r2)
    {
        nbh.push_back(pointIndex(*branch));
    }

    //
    KDNodePtr section;
    KDNodePtr other;
    if (dx > 0)
    {
        section = branch->m_left;
        other = branch->m_right;
    }
    else
    {
        section = branch->m_right;
        other = branch->m_left;
    }

    nbh_s = neighborhood_(section, pt, rad, (level + 1) % dim);
    nbh.insert(nbh.end(), nbh_s.begin(), nbh_s.end());
    if (dx2 < r2)
    {
        nbh_o = neighborhood_(other, pt, rad, (level + 1) % dim);
        nbh.insert(nbh.end(), nbh_o.begin(), nbh_o.end());
    }

    return nbh;
}

point_t KDTree::nearest_point(const point_t &pt)
{
    return point_t(*nearest_(pt));
}

size_t KDTree::nearest_index(const point_t &pt)
{
    return size_t(*nearest_(pt));
}

pointIndex KDTree::nearest_pointIndex(const point_t &pt)
{
    KDNodePtr Nearest = nearest_(pt);
    return pointIndex(point_t(*Nearest), size_t(*Nearest));
}

pointIndexArr KDTree::neighborhood(const point_t& pt, const double& rad)
{
    size_t level = 0;
    return neighborhood_(m_root, pt, rad, level);
}

pointVec KDTree::neighborhood_points(const point_t& pt, const double& rad)
{
    size_t level = 0;
    pointIndexArr nbh = neighborhood_(m_root, pt, rad, level);
    pointVec nbhp;
    nbhp.resize(nbh.size());
    std::transform(nbh.begin(), nbh.end(), nbhp.begin(),
                   [](pointIndex x) { return x.first; });
    return nbhp;
}

indexArr KDTree::neighborhood_indices(const point_t& pt, const double& rad)
{
    size_t level = 0;
    pointIndexArr nbh = neighborhood_(m_root, pt, rad, level);
    indexArr nbhi;
    nbhi.resize(nbh.size());
    std::transform(nbh.begin(), nbh.end(), nbhi.begin(),
                   [](pointIndex x) { return x.second; });
    return nbhi;
}
