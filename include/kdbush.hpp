// https://github.com/mourner/kdbush.hpp

#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

template<class TFloat = double, class TIndex = std::size_t>
class KDBush {
    
public:

    static const uint8_t defaultNodeSize = 64;
    
    KDBush(const std::vector<TFloat>& coords, const uint8_t nodeSize_ = defaultNodeSize)
        : KDBush(std::begin(coords), std::end(coords), nodeSize_) {}

    template<class TFloatIter>
    KDBush(TFloatIter coordsBegin, TFloatIter coordsEnd, const uint8_t nodeSize_ = defaultNodeSize)
        : coords(std::distance(coordsBegin, coordsEnd)), nodeSize(nodeSize_) {

        std::copy(coordsBegin, coordsEnd, coords.data());

        const auto ids_size = coords.size() / 2;
        ids.reserve(ids_size);
        for (size_t i = 0; i < ids_size; i++) ids.push_back(TIndex(i));

        sortKD(0, ids.size() - 1, 0);
    }

    void range(const TFloat minX,
               const TFloat minY,
               const TFloat maxX,
               const TFloat maxY,
               std::vector<TIndex> &result) {
        range(minX, minY, maxX, maxY, [&result] (TIndex index) { result.push_back(index); });
    }

    template<class TPushBackIndex>
    void range(const TFloat minX,
               const TFloat minY,
               const TFloat maxX,
               const TFloat maxY,
               const TPushBackIndex& push_back_index) {
        range(minX, minY, maxX, maxY, push_back_index, 0, ids.size() - 1, 0);
    }

    void within(const TFloat qx, const TFloat qy, const TFloat r, std::vector<TIndex>& result) {
        within(qx, qy, r, [&result] (TIndex index) { result.push_back(index); });
    }

    template<class TPushBackIndex>
    void within(const TFloat qx, const TFloat qy, const TFloat r, const TPushBackIndex& push_back_index) {
        within(qx, qy, r, push_back_index, 0, ids.size() - 1, 0);
    }

private:
    uint8_t nodeSize;
    std::vector<TIndex> ids;
    std::vector<TFloat> coords;

    template<class TPushBackIndex>
    void range(const TFloat minX,
               const TFloat minY,
               const TFloat maxX,
               const TFloat maxY,
               const TPushBackIndex& push_back_index,
               const TIndex left,
               const TIndex right,
               const uint8_t axis) {

        if (right - left <= nodeSize) {
            for (auto i = left; i <= right; i++) {
                const TFloat x = coords[2 * i];
                const TFloat y = coords[2 * i + 1];
                if (x >= minX && x <= maxX && y >= minY && y <= maxY) push_back_index(ids[i]);
            }
            return;
        }

        const TIndex m = (left + right) >> 1;
        const TFloat x = coords[2 * m];
        const TFloat y = coords[2 * m + 1];

        if (x >= minX && x <= maxX && y >= minY && y <= maxY) push_back_index(ids[m]);

        if (axis == 0 ? minX <= x : minY <= y)
            range(minX, minY, maxX, maxY, result, left, m - 1, (axis + 1) % 2);

        if (axis == 0 ? maxX >= x : maxY >= y)
            range(minX, minY, maxX, maxY, result, m + 1, right, (axis + 1) % 2);
    }

    template<class TPushBackIndex>    
    void within(const TFloat qx,
                const TFloat qy,
                const TFloat r,
                const TPushBackIndex& push_back_index,
                const TIndex left,
                const TIndex right,
                const uint8_t axis) {

        const TFloat r2 = r * r;

        if (right - left <= nodeSize) {
            for (auto i = left; i <= right; i++) {
                if (sqDist(coords[2 * i], coords[2 * i + 1], qx, qy) <= r2)
                    push_back_index(ids[i]);
            }
            return;
        }

        const TIndex m = (left + right) >> 1;
        const TFloat x = coords[2 * m];
        const TFloat y = coords[2 * m + 1];

        if (sqDist(x, y, qx, qy) <= r2) push_back_index(ids[m]);

        if (axis == 0 ? qx - r <= x : qy - r <= y)
            within(qx, qy, r, push_back_index, left, m - 1, (axis + 1) % 2);

        if (axis == 0 ? qx + r >= x : qy + r >= y)
            within(qx, qy, r, push_back_index, m + 1, right, (axis + 1) % 2);
    }

    void sortKD(const TIndex left, const TIndex right, const uint8_t axis) {
        if (right - left <= nodeSize) return;
        const TIndex m = (left + right) >> 1;
        select(m, left, right, axis);
        sortKD(left, m - 1, (axis + 1) % 2);
        sortKD(m + 1, right, (axis + 1) % 2);
    }

    void select(const TIndex k, TIndex left, TIndex right, const uint8_t axis) {

        while (right > left) {
            if (right - left > 600) {
                const TIndex n = right - left + 1;
                const TIndex m = k - left + 1;
                const TFloat z = log(n);
                const TFloat s = 0.5 * exp(2 * z / 3);
                const TFloat sd = 0.5 * sqrt(z * s * (n - s) / n) * (2 * m < n ? -1 : 1);
                const TIndex newLeft = std::max(left, TIndex(k - m * s / n + sd));
                const TIndex newRight = std::min(right, TIndex(k + (n - m) * s / n + sd));
                select(k, newLeft, newRight, axis);
            }

            const TFloat t = coords[2 * k + axis];
            TIndex i = left;
            TIndex j = right;

            swapItem(left, k);
            if (coords[2 * right + axis] > t) swapItem(left, right);

            while (i < j) {
                swapItem(i, j);
                i++;
                j--;
                while (coords[2 * i + axis] < t) i++;
                while (coords[2 * j + axis] > t) j--;
            }

            if (coords[2 * left + axis] == t)
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
        std::iter_swap(coords.begin() + 2 * i, coords.begin() + 2 * j);
        std::iter_swap(coords.begin() + 2 * i + 1, coords.begin() + 2 * j + 1);
    }

    TFloat sqDist(const TFloat ax, const TFloat ay, const TFloat bx, const TFloat by) {
        const TFloat dx = ax - bx;
        const TFloat dy = ay - by;
        return dx * dx + dy * dy;
    }
};
