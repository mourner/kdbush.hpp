#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <tuple>
#include <vector>

namespace kdbush {

template <std::uint8_t I, typename T>
struct nth {
    inline static typename std::tuple_element<I, T>::type get(const T &t) {
        return std::get<I>(t);
    }
};

template <typename TPoint, typename TIndex = std::size_t>
class KDBush {

public:
    using TNumber = decltype(nth<0, TPoint>::get(std::declval<TPoint>()));
    static_assert(
        std::is_same<TNumber, decltype(nth<1, TPoint>::get(std::declval<TPoint>()))>::value,
        "point component types must be identical");

    static const std::uint8_t defaultNodeSize = 64;

    KDBush(const std::vector<TPoint> &points_, const std::uint8_t nodeSize_ = defaultNodeSize)
        : KDBush(std::begin(points_), std::end(points_), nodeSize_) {
    }

    template <typename TPointIter>
    KDBush(const TPointIter &points_begin,
           const TPointIter &points_end,
           const std::uint8_t nodeSize_ = defaultNodeSize)
        : nodeSize(nodeSize_) {

        const TIndex size = std::distance(points_begin, points_end);

        points.reserve(size);
        ids.reserve(size);

        for (TIndex i = 0; i < size; i++) {
            const auto p = *(points_begin + i);
            points.emplace_back(nth<0, TPoint>::get(p), nth<1, TPoint>::get(p));
            ids.push_back(i);
        }

        sortKD(0, size - 1, 0);
    }

    template <typename TOutputIter>
    void range(const TNumber minX,
               const TNumber minY,
               const TNumber maxX,
               const TNumber maxY,
               TOutputIter out) {
        range(minX, minY, maxX, maxY, out, 0, ids.size() - 1, 0);
    }

    template <typename TOutputIter>
    void within(const TNumber qx, const TNumber qy, const TNumber r, TOutputIter out) {
        within(qx, qy, r, out, 0, ids.size() - 1, 0);
    }

private:
    std::vector<TIndex> ids;
    std::vector<std::pair<TNumber, TNumber>> points;
    std::uint8_t nodeSize;

    template <typename TOutputIter>
    void range(const TNumber minX,
               const TNumber minY,
               const TNumber maxX,
               const TNumber maxY,
               TOutputIter out,
               const TIndex left,
               const TIndex right,
               const std::uint8_t axis) {

        if (right - left <= nodeSize) {
            for (auto i = left; i <= right; i++) {
                const TNumber x = std::get<0>(points[i]);
                const TNumber y = std::get<1>(points[i]);
                if (x >= minX && x <= maxX && y >= minY && y <= maxY) *out++ = ids[i];
            }
            return;
        }

        const TIndex m = (left + right) >> 1;
        const TNumber x = std::get<0>(points[m]);
        const TNumber y = std::get<1>(points[m]);

        if (x >= minX && x <= maxX && y >= minY && y <= maxY) *out++ = ids[m];

        if (axis == 0 ? minX <= x : minY <= y)
            range(minX, minY, maxX, maxY, out, left, m - 1, (axis + 1) % 2);

        if (axis == 0 ? maxX >= x : maxY >= y)
            range(minX, minY, maxX, maxY, out, m + 1, right, (axis + 1) % 2);
    }

    template <typename TOutputIter>
    void within(const TNumber qx,
                const TNumber qy,
                const TNumber r,
                TOutputIter out,
                const TIndex left,
                const TIndex right,
                const std::uint8_t axis) {

        const TNumber r2 = r * r;

        if (right - left <= nodeSize) {
            for (auto i = left; i <= right; i++) {
                const TNumber x = std::get<0>(points[i]);
                const TNumber y = std::get<1>(points[i]);
                if (sqDist(x, y, qx, qy) <= r2) *out++ = ids[i];
            }
            return;
        }

        const TIndex m = (left + right) >> 1;
        const TNumber x = std::get<0>(points[m]);
        const TNumber y = std::get<1>(points[m]);

        if (sqDist(x, y, qx, qy) <= r2) *out++ = ids[m];

        if (axis == 0 ? qx - r <= x : qy - r <= y)
            within(qx, qy, r, out, left, m - 1, (axis + 1) % 2);

        if (axis == 0 ? qx + r >= x : qy + r >= y)
            within(qx, qy, r, out, m + 1, right, (axis + 1) % 2);
    }

    void sortKD(const TIndex left, const TIndex right, const std::uint8_t axis) {
        if (right - left <= nodeSize) return;
        const TIndex m = (left + right) >> 1;
        if (axis == 0) {
            select<0>(m, left, right);
        } else {
            select<1>(m, left, right);
        }
        sortKD(left, m - 1, (axis + 1) % 2);
        sortKD(m + 1, right, (axis + 1) % 2);
    }

    template <std::uint8_t axis>
    void select(const TIndex k, TIndex left, TIndex right) {

        while (right > left) {
            if (right - left > 600) {
                const TIndex n = right - left + 1;
                const TIndex m = k - left + 1;
                const double z = log(n);
                const double s = 0.5 * exp(2 * z / 3);
                const double sd = 0.5 * sqrt(z * s * (n - s) / n) * (2 * m < n ? -1 : 1);
                const TIndex newLeft = std::max(left, TIndex(k - m * s / n + sd));
                const TIndex newRight = std::min(right, TIndex(k + (n - m) * s / n + sd));
                select<axis>(k, newLeft, newRight);
            }

            const TNumber t = std::get<axis>(points[k]);
            TIndex i = left;
            TIndex j = right;

            swapItem(left, k);
            if (std::get<axis>(points[right]) > t) swapItem(left, right);

            while (i < j) {
                swapItem(i, j);
                i++;
                j--;
                while (std::get<axis>(points[i]) < t) i++;
                while (std::get<axis>(points[j]) > t) j--;
            }

            if (std::get<axis>(points[left]) == t)
                swapItem(left, j);
            else {
                j++;
                swapItem(j, right);
            }

            if (j <= k) left = j + 1;
            if (k <= j) right = j - 1;
        }
    }

    void swapItem(const TIndex i, const TIndex j) {
        std::iter_swap(ids.begin() + i, ids.begin() + j);
        std::iter_swap(points.begin() + i, points.begin() + j);
    }

    TNumber sqDist(const TNumber ax, const TNumber ay, const TNumber bx, const TNumber by) {
        return std::pow(ax - bx, 2) + std::pow(ay - by, 2);
    }
};

} // namespace kdbush
