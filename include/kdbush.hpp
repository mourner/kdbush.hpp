#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

class KDBush {

    using size_t = std::size_t;

public:
    KDBush(const std::vector<double> &coords_, const uint8_t nodeSize_ = 64)
        : coords(coords_), nodeSize(nodeSize_) {

        const auto ids_size = coords.size() / 2;
        ids.reserve(ids_size);
        for (size_t i = 0; i < ids_size; i++) ids.push_back(i);

        sortKD(0, ids.size() - 1, 0);
    }

    void within(const double qx,
                const double qy,
                const double r,
                std::vector<size_t> &result,
                const size_t left,
                const size_t right,
                const uint8_t axis) {

        const double r2 = r * r;

        if (right - left <= nodeSize) {
            for (auto i = left; i <= right; i++) {
                if (sqDist(coords[2 * i], coords[2 * i + 1], qx, qy) <= r2)
                    result.push_back(ids[i]);
            }
            return;
        }

        const size_t m = (left + right) >> 1;
        const double x = coords[2 * m];
        const double y = coords[2 * m + 1];

        if (sqDist(x, y, qx, qy) <= r2) result.push_back(ids[m]);

        if (axis == 0 ? qx - r <= x : qy - r <= y)
            within(qx, qy, r, result, left, m - 1, (axis + 1) % 2);

        if (axis == 0 ? qx + r >= x : qy + r >= y)
            within(qx, qy, r, result, m + 1, right, (axis + 1) % 2);
    }

    void within(const double qx, const double qy, const double r, std::vector<size_t> &result) {
        within(qx, qy, r, result, 0, ids.size() - 1, 0);
    }

private:
    uint8_t nodeSize;
    std::vector<size_t> ids;
    std::vector<double> coords;

    void sortKD(const size_t left, const size_t right, const uint8_t axis) {
        if (right - left <= nodeSize) return;
        const size_t m = (left + right) >> 1;
        select(m, left, right, axis);
        sortKD(left, m - 1, (axis + 1) % 2);
        sortKD(m + 1, right, (axis + 1) % 2);
    }

    void select(const size_t k, size_t left, size_t right, const uint8_t axis) {

        while (right > left) {
            if (right - left > 600) {
                const size_t n = right - left + 1;
                const size_t m = k - left + 1;
                const double z = log(n);
                const double s = 0.5 * exp(2 * z / 3);
                const double sd = 0.5 * sqrt(z * s * (n - s) / n) * (2 * m < n ? -1 : 1);
                const size_t newLeft = std::max(left, size_t(k - m * s / n + sd));
                const size_t newRight = std::min(right, size_t(k + (n - m) * s / n + sd));
                select(k, newLeft, newRight, axis);
            }

            const double t = coords[2 * k + axis];
            size_t i = left;
            size_t j = right;

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

    void swapItem(const size_t i, const size_t j) {
        std::iter_swap(ids.begin() + i, ids.begin() + j);
        std::iter_swap(coords.begin() + 2 * i, coords.begin() + 2 * j);
        std::iter_swap(coords.begin() + 2 * i + 1, coords.begin() + 2 * j + 1);
    }

    double sqDist(const double ax, const double ay, const double bx, const double by) {
        const double dx = ax - bx;
        const double dy = ay - by;
        return dx * dx + dy * dy;
    }
};
