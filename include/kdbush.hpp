#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

template<class TNumber = double, class TIndex = std::size_t>
class KDBush {

public:
    KDBush(const std::vector<TNumber> &coords_, const uint8_t nodeSize_ = 64)
        : coords(coords_), nodeSize(nodeSize_) {

        const auto ids_size = coords.size() / 2;
        ids.reserve(ids_size);
        for (TIndex i = 0; i < ids_size; i++) ids.push_back(i);

        sortKD(0, ids.size() - 1, 0);
    }

    void range(const TNumber minX,
               const TNumber minY,
               const TNumber maxX,
               const TNumber maxY,
               std::vector<TIndex> &result) {
        range(minX, minY, maxX, maxY, result, 0, ids.size() - 1, 0);
    }

    void within(const TNumber qx, const TNumber qy, const TNumber r, std::vector<TIndex> &result) {
        within(qx, qy, r, result, 0, ids.size() - 1, 0);
    }

private:
    std::vector<TIndex> ids;
    std::vector<TNumber> coords;
    uint8_t nodeSize;

    void range(const TNumber minX,
               const TNumber minY,
               const TNumber maxX,
               const TNumber maxY,
               std::vector<TIndex> &result,
               const TIndex left,
               const TIndex right,
               const uint8_t axis) {

        if (right - left <= nodeSize) {
            for (auto i = left; i <= right; i++) {
                const TNumber x = coords[2 * i];
                const TNumber y = coords[2 * i + 1];
                if (x >= minX && x <= maxX && y >= minY && y <= maxY) result.push_back(ids[i]);
            }
            return;
        }

        const TIndex m = (left + right) >> 1;
        const TNumber x = coords[2 * m];
        const TNumber y = coords[2 * m + 1];

        if (x >= minX && x <= maxX && y >= minY && y <= maxY) result.push_back(ids[m]);

        if (axis == 0 ? minX <= x : minY <= y)
            range(minX, minY, maxX, maxY, result, left, m - 1, (axis + 1) % 2);

        if (axis == 0 ? maxX >= x : maxY >= y)
            range(minX, minY, maxX, maxY, result, m + 1, right, (axis + 1) % 2);
    }

    void within(const TNumber qx,
                const TNumber qy,
                const TNumber r,
                std::vector<TIndex> &result,
                const TIndex left,
                const TIndex right,
                const uint8_t axis) {

        const TNumber r2 = r * r;

        if (right - left <= nodeSize) {
            for (auto i = left; i <= right; i++) {
                if (sqDist(coords[2 * i], coords[2 * i + 1], qx, qy) <= r2)
                    result.push_back(ids[i]);
            }
            return;
        }

        const TIndex m = (left + right) >> 1;
        const TNumber x = coords[2 * m];
        const TNumber y = coords[2 * m + 1];

        if (sqDist(x, y, qx, qy) <= r2) result.push_back(ids[m]);

        if (axis == 0 ? qx - r <= x : qy - r <= y)
            within(qx, qy, r, result, left, m - 1, (axis + 1) % 2);

        if (axis == 0 ? qx + r >= x : qy + r >= y)
            within(qx, qy, r, result, m + 1, right, (axis + 1) % 2);
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
                const double z = log(n);
                const double s = 0.5 * exp(2 * z / 3);
                const double sd = 0.5 * sqrt(z * s * (n - s) / n) * (2 * m < n ? -1 : 1);
                const TIndex newLeft = std::max(left, TIndex(k - m * s / n + sd));
                const TIndex newRight = std::min(right, TIndex(k + (n - m) * s / n + sd));
                select(k, newLeft, newRight, axis);
            }

            const TNumber t = coords[2 * k + axis];
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

    TNumber sqDist(const TNumber ax, const TNumber ay, const TNumber bx, const TNumber by) {
        const TNumber dx = ax - bx;
        const TNumber dy = ay - by;
        return dx * dx + dy * dy;
    }
};
