#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

class KDBush {
    uint32_t nodeSize;
    std::vector<double> coords;

    void sortKD(uint32_t left, uint32_t right, uint8_t axis) {
        if (right - left <= nodeSize) return;
        uint32_t m = (left + right) >> 1;
        select(m, left, right, axis);
        sortKD(left, m - 1, (axis + 1) % 2);
        sortKD(m + 1, right, (axis + 1) % 2);
    }

    void select(uint32_t k, uint32_t left, uint32_t right, uint8_t axis) {

        while (right > left) {
            if (right - left > 600) {
                uint32_t n = right - left + 1;
                uint32_t m = k - left + 1;
                double z = log(n);
                double s = 0.5 * exp(2 * z / 3);
                double sd = 0.5 * sqrt(z * s * (n - s) / n);
                if (2 * m < n) sd = -sd;
                uint32_t newLeft = std::max(left, uint32_t(k - m * s / n + sd));
                uint32_t newRight = std::min(right, uint32_t(k + (n - m) * s / n + sd));
                select(k, newLeft, newRight, axis);
            }

            auto t = coords[2 * k + axis];
            uint32_t i = left;
            uint32_t j = right;

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

    void swapItem(uint32_t i, uint32_t j) {
        std::iter_swap(ids.begin() + i, ids.begin() + j);
        std::iter_swap(coords.begin() + 2 * i, coords.begin() + 2 * j);
        std::iter_swap(coords.begin() + 2 * i + 1, coords.begin() + 2 * j + 1);
    }

    void search_within(uint32_t left,
                       uint32_t right,
                       uint8_t axis,
                       double qx,
                       double qy,
                       double r,
                       std::vector<uint32_t> &result) {

        double r2 = r * r;

        if (right - left <= nodeSize) {
            for (auto i = left; i <= right; i++) {
                if (sqDist(coords[2 * i], coords[2 * i + 1], qx, qy) <= r2)
                    result.push_back(ids[i]);
            }
            return;
        }

        uint32_t m = (left + right) >> 1;
        double x = coords[2 * m];
        double y = coords[2 * m + 1];

        if (sqDist(x, y, qx, qy) <= r2) result.push_back(ids[m]);

        if (axis == 0 ? qx - r <= x : qy - r <= y)
            search_within(left, m - 1, (axis + 1) % 2, qx, qy, r, result);

        if (axis == 0 ? qx + r >= x : qy + r >= y)
            search_within(m + 1, right, (axis + 1) % 2, qx, qy, r, result);
    }

    double sqDist(double ax, double ay, double bx, double by) {
        double dx = ax - bx;
        double dy = ay - by;
        return dx * dx + dy * dy;
    }

public:
    std::vector<uint32_t> ids;

    KDBush(std::vector<double> const &coords_, uint8_t nodeSize_ = 64)
        : coords(coords_), nodeSize(nodeSize_) {
        uint32_t ids_size = coords.size() / 2;
        ids.reserve(ids_size);
        for (uint32_t i = 0; i < ids_size; i++) ids.push_back(i);

        sortKD(0, ids.size() - 1, 0);
    }

    std::vector<uint32_t> within(double qx, double qy, double r) {
        std::vector<uint32_t> result;
        search_within(0, ids.size() - 1, 0, qx, qy, r, result);
        return result;
    }
};
