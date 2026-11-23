/**************************************************************************/
/*  delaunay_3d.h                                                         */
/**************************************************************************/
/*                         This file is part of:                          */
/*                             GODOT ENGINE                               */
/*                        https://godotengine.org                         */
/**************************************************************************/
/* Copyright (c) 2014-present Godot Engine contributors (see AUTHORS.md). */
/* Copyright (c) 2007-2014 Juan Linietsky, Ariel Manzur.                  */
/*                                                                        */
/* Permission is hereby granted, free of charge, to any person obtaining  */
/* a copy of this software and associated documentation files (the        */
/* "Software"), to deal in the Software without restriction, including    */
/* without limitation the rights to use, copy, modify, merge, publish,    */
/* distribute, sublicense, and/or sell copies of the Software, and to     */
/* permit persons to whom the Software is furnished to do so, subject to  */
/* the following conditions:                                              */
/*                                                                        */
/* The above copyright notice and this permission notice shall be         */
/* included in all copies or substantial portions of the Software.        */
/*                                                                        */
/* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,        */
/* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF     */
/* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. */
/* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY   */
/* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,   */
/* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE      */
/* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                 */
/**************************************************************************/

#pragma once

#include "core/math/aabb.h"
#include "core/math/projection.h"
#include "core/math/vector3.h"
#include "core/math/vector3i.h"
#include "core/templates/a_hash_map.h"
#include "core/templates/list.h"
#include "core/templates/local_vector.h"
#include "core/templates/vector.h"

class Delaunay3D {
	struct Simplex;

	enum {
		ACCEL_GRID_SIZE = 16,
		QUANTIZATION_MAX = 1 << 16 // A power of two smaller than the 23 bit significand of a float.
	};
	struct GridPos {
		Vector3i pos;
		List<Simplex *>::Element *E = nullptr;
	};

	struct Simplex {
		uint32_t points[4];
		int orientation = 0; // Sign of the exact orientation of the four points.
		LocalVector<GridPos> grid_positions;
		List<Simplex *>::Element *SE = nullptr;

		_FORCE_INLINE_ Simplex() {}
		_FORCE_INLINE_ Simplex(uint32_t p_a, uint32_t p_b, uint32_t p_c, uint32_t p_d) {
			points[0] = p_a;
			points[1] = p_b;
			points[2] = p_c;
			points[3] = p_d;
		}
	};

	struct Triangle {
		uint32_t triangle[3];
		bool bad = false;
		_FORCE_INLINE_ bool operator==(const Triangle &p_triangle) const {
			return triangle[0] == p_triangle.triangle[0] && triangle[1] == p_triangle.triangle[1] && triangle[2] == p_triangle.triangle[2];
		}

		_FORCE_INLINE_ Triangle() {}
		_FORCE_INLINE_ Triangle(uint32_t p_a, uint32_t p_b, uint32_t p_c) {
			if (p_a > p_b) {
				SWAP(p_a, p_b);
			}
			if (p_b > p_c) {
				SWAP(p_b, p_c);
			}
			if (p_a > p_b) {
				SWAP(p_a, p_b);
			}

			triangle[0] = p_a;
			triangle[1] = p_b;
			triangle[2] = p_c;
		}
	};

	struct TriangleHasher {
		_FORCE_INLINE_ static uint32_t hash(const Triangle &p_triangle) {
			uint32_t h = hash_djb2_one_32(p_triangle.triangle[0]);
			h = hash_djb2_one_32(p_triangle.triangle[1], h);
			return hash_fmix32(hash_djb2_one_32(p_triangle.triangle[2], h));
		}
	};

	// The predicates below work on the quantized integer coordinates, so they can be evaluated exactly.
	// Points are quantized to [0, QUANTIZATION_MAX] and the root simplex is about 141 times larger,
	// so coordinate differences stay below 2^25 and the in-sphere determinant below 2^132.
	// A double evaluation with an error bound decides almost all cases, the rest is computed exactly
	// with a small fixed size 256-bit integer.

	static constexpr double UNIT_ROUNDOFF = 1.0 / 9007199254740992.0; // 2^-53.

	// Signed 256-bit two's complement integer, only with the operations the predicates need.
	struct WideInt {
		uint64_t w[4] = { 0, 0, 0, 0 };

		_FORCE_INLINE_ WideInt() {}
		_FORCE_INLINE_ explicit WideInt(int64_t p_value) {
			const uint64_t fill = p_value < 0 ? ~uint64_t(0) : 0;
			w[0] = uint64_t(p_value);
			w[1] = fill;
			w[2] = fill;
			w[3] = fill;
		}

		_FORCE_INLINE_ WideInt operator+(const WideInt &p_other) const {
			WideInt r;
			uint64_t carry = 0;
			for (int i = 0; i < 4; i++) {
				const uint64_t a = w[i] + carry;
				const uint64_t carry_a = a < carry;
				r.w[i] = a + p_other.w[i];
				carry = carry_a + (r.w[i] < a);
			}
			return r;
		}

		_FORCE_INLINE_ WideInt operator-() const {
			WideInt r;
			for (int i = 0; i < 4; i++) {
				r.w[i] = ~w[i];
			}
			return r + WideInt(1);
		}

		_FORCE_INLINE_ WideInt operator-(const WideInt &p_other) const {
			return *this + (-p_other);
		}

		_FORCE_INLINE_ bool is_negative() const {
			return w[3] >> 63;
		}

		_FORCE_INLINE_ int sign() const {
			if (is_negative()) {
				return -1;
			}
			return (w[0] | w[1] | w[2] | w[3]) ? 1 : 0;
		}

		// Full 64x64 -> 128 bit unsigned product, portable.
		_FORCE_INLINE_ static void mul_64(uint64_t p_a, uint64_t p_b, uint64_t &r_lo, uint64_t &r_hi) {
			const uint64_t a0 = p_a & 0xFFFFFFFF;
			const uint64_t a1 = p_a >> 32;
			const uint64_t b0 = p_b & 0xFFFFFFFF;
			const uint64_t b1 = p_b >> 32;
			const uint64_t p00 = a0 * b0;
			const uint64_t p01 = a0 * b1;
			const uint64_t p10 = a1 * b0;
			const uint64_t p11 = a1 * b1;
			const uint64_t mid = (p00 >> 32) + (p01 & 0xFFFFFFFF) + (p10 & 0xFFFFFFFF);
			r_lo = (p00 & 0xFFFFFFFF) | (mid << 32);
			r_hi = p11 + (p01 >> 32) + (p10 >> 32) + (mid >> 32);
		}

		// The magnitudes used here stay far below 2^255, so the result never overflows.
		_FORCE_INLINE_ WideInt operator*(int64_t p_b) const {
			const bool negative = is_negative() != (p_b < 0);
			const WideInt a = is_negative() ? -*this : *this;
			const uint64_t b = p_b < 0 ? uint64_t(0) - uint64_t(p_b) : uint64_t(p_b);
			WideInt r;
			uint64_t carry = 0;
			for (int i = 0; i < 4; i++) {
				uint64_t lo, hi;
				mul_64(a.w[i], b, lo, hi);
				lo += carry;
				hi += lo < carry;
				r.w[i] = lo;
				carry = hi;
			}
			return negative ? -r : r;
		}
	};

	_FORCE_INLINE_ static WideInt det3_exact(const int64_t *p_a, const int64_t *p_b, const int64_t *p_c) {
		return (WideInt(p_b[1]) * p_c[2] - WideInt(p_b[2]) * p_c[1]) * p_a[0] - (WideInt(p_b[0]) * p_c[2] - WideInt(p_b[2]) * p_c[0]) * p_a[1] + (WideInt(p_b[0]) * p_c[1] - WideInt(p_b[1]) * p_c[0]) * p_a[2];
	}

	_FORCE_INLINE_ static double det3(const double *p_a, const double *p_b, const double *p_c) {
		return p_a[0] * (p_b[1] * p_c[2] - p_b[2] * p_c[1]) - p_a[1] * (p_b[0] * p_c[2] - p_b[2] * p_c[0]) + p_a[2] * (p_b[0] * p_c[1] - p_b[1] * p_c[0]);
	}

	// Same products as det3(), but every term added with its absolute value. Used for the error bound.
	_FORCE_INLINE_ static double det3_permanent(const double *p_a, const double *p_b, const double *p_c) {
		return p_a[0] * (p_b[1] * p_c[2] + p_b[2] * p_c[1]) + p_a[1] * (p_b[0] * p_c[2] + p_b[2] * p_c[0]) + p_a[2] * (p_b[0] * p_c[1] + p_b[1] * p_c[0]);
	}

	// Sign of det(b - a, c - a, d - a).
	static int orientation_sign(const Vector3i *p_points, const Simplex &p_simplex) {
		const Vector3i &a = p_points[p_simplex.points[0]];
		int64_t rows[3][3];
		for (int i = 0; i < 3; i++) {
			const Vector3i &p = p_points[p_simplex.points[i + 1]];
			for (int k = 0; k < 3; k++) {
				rows[i][k] = int64_t(p[k]) - int64_t(a[k]);
			}
		}
		return det3_exact(rows[0], rows[1], rows[2]).sign();
	}

	// Sign of the in-sphere determinant of the simplex relative to p_vertex.
	// For a positively oriented simplex, a positive value means p_vertex is strictly inside its circumsphere.
	static int insphere_sign(const Vector3i *p_points, const Simplex &p_simplex, uint32_t p_vertex) {
		const Vector3i &e = p_points[p_vertex];
		int64_t rows[4][4];
		double rows_d[4][4];
		double rows_abs[4][4];
		for (int i = 0; i < 4; i++) {
			const Vector3i &p = p_points[p_simplex.points[i]];
			for (int k = 0; k < 3; k++) {
				rows[i][k] = int64_t(p[k]) - int64_t(e[k]);
				rows_d[i][k] = double(rows[i][k]);
				rows_abs[i][k] = Math::abs(rows_d[i][k]);
			}
			// Exact, since the sum of the squares stays below 2^53.
			rows_d[i][3] = rows_d[i][0] * rows_d[i][0] + rows_d[i][1] * rows_d[i][1] + rows_d[i][2] * rows_d[i][2];
			rows_abs[i][3] = rows_d[i][3];
			rows[i][3] = int64_t(rows_d[i][3]);
		}

		// Fast path in double, with a strict error bound. Let u = 2^-53 be the unit roundoff.
		// The coordinate differences are integers below 2^25, so the 2x2 minors inside det3() and the
		// lifted coordinates are integers below 2^53 and are computed exactly. Rounding only happens
		// after that: one product and two sums in det3(), one product and three sums outside. Each term
		// goes through at most 7 roundings, so the error is at most 7u / (1 - 7u) times the permanent
		// (the same expression on absolute values, see det3_permanent()). The computed permanent has
		// the same relative error, so the error is below 7.01u times the computed permanent, and
		// 8u = 2^-50 covers it. Multiplying by a power of two is exact. Fused multiply-add only
		// removes roundings, so the bound also holds if the compiler contracts the expressions.
		const double value = rows_d[0][3] * det3(rows_d[1], rows_d[2], rows_d[3]) - rows_d[1][3] * det3(rows_d[0], rows_d[2], rows_d[3]) + rows_d[2][3] * det3(rows_d[0], rows_d[1], rows_d[3]) - rows_d[3][3] * det3(rows_d[0], rows_d[1], rows_d[2]);
		const double permanent = rows_abs[0][3] * det3_permanent(rows_abs[1], rows_abs[2], rows_abs[3]) + rows_abs[1][3] * det3_permanent(rows_abs[0], rows_abs[2], rows_abs[3]) + rows_abs[2][3] * det3_permanent(rows_abs[0], rows_abs[1], rows_abs[3]) + rows_abs[3][3] * det3_permanent(rows_abs[0], rows_abs[1], rows_abs[2]);
		const double bound = permanent * (8.0 * UNIT_ROUNDOFF);
		if (value > bound) {
			return 1;
		}
		if (value < -bound) {
			return -1;
		}

		// Too close to call, compute exactly.
		const WideInt exact = det3_exact(rows[1], rows[2], rows[3]) * rows[0][3] - det3_exact(rows[0], rows[2], rows[3]) * rows[1][3] + det3_exact(rows[0], rows[1], rows[3]) * rows[2][3] - det3_exact(rows[0], rows[1], rows[2]) * rows[3][3];
		return exact.sign();
	}

	// Bounding sphere of the circumsphere, in quantized coordinates. It is used to register the simplex
	// in the acceleration grid, so it must contain the whole circumsphere: if a simplex whose circumsphere
	// contains a new point is not found, the cavity is wrong and the triangulation breaks.
	// Returns false if the bound can't be established (almost flat simplex), then the whole grid is used.
	static bool circum_sphere_bound(const Vector3i *p_points, const Simplex &p_simplex, double *r_center, double &r_radius) {
		const Vector3i &origin = p_points[p_simplex.points[0]];
		double cols[3][3]; // Edge vectors from the origin vertex, as columns.
		double cols_abs[3][3];
		double sq[3]; // Squared edge lengths.
		for (int i = 0; i < 3; i++) {
			const Vector3i &p = p_points[p_simplex.points[i + 1]];
			for (int k = 0; k < 3; k++) {
				cols[k][i] = double(int64_t(p[k]) - int64_t(origin[k]));
				cols_abs[k][i] = Math::abs(cols[k][i]);
			}
		}
		for (int i = 0; i < 3; i++) {
			sq[i] = cols[0][i] * cols[0][i] + cols[1][i] * cols[1][i] + cols[2][i] * cols[2][i];
		}

		// The center relative to the origin vertex is N / (2 D) by Cramer's rule, where D is the
		// determinant of the edge vectors and N[k] is the same determinant with column k replaced by sq.
		// The row passed first to det3() is the only one not used in its 2x2 minors, so sq always goes
		// first (cyclic permutations keep the sign). Then, like in insphere_sign(), the minors are exact
		// and each term goes through at most 3 roundings: the error is below 3.01u times the computed
		// permanent, 4u = 2^-51 covers it.
		const double d = det3(cols[0], cols[1], cols[2]);
		const double d_error = det3_permanent(cols_abs[0], cols_abs[1], cols_abs[2]) * (4.0 * UNIT_ROUNDOFF);
		if (!(Math::abs(d) > 2.0 * d_error)) {
			return false;
		}
		const double n[3] = {
			det3(sq, cols[1], cols[2]),
			det3(sq, cols[2], cols[0]),
			det3(sq, cols[0], cols[1]),
		};
		const double n_error[3] = {
			det3_permanent(sq, cols_abs[1], cols_abs[2]) * (4.0 * UNIT_ROUNDOFF),
			det3_permanent(sq, cols_abs[2], cols_abs[0]) * (4.0 * UNIT_ROUNDOFF),
			det3_permanent(sq, cols_abs[0], cols_abs[1]) * (4.0 * UNIT_ROUNDOFF),
		};

		// For the true values: |N/D - n/d| <= (|N - n| + |n/d| |D - d|) / (|d| - |D - d|).
		// The division adds one more rounding of the computed center.
		double center_error = 0.0;
		double center_length2 = 0.0;
		for (int k = 0; k < 3; k++) {
			const double c = n[k] / (2.0 * d);
			const double error = (n_error[k] + 2.0 * Math::abs(c) * d_error) / (2.0 * (Math::abs(d) - d_error)) + 2.0 * UNIT_ROUNDOFF * Math::abs(c);
			center_error += error;
			center_length2 += c * c;
			r_center[k] = c;
		}

		// The true sphere is inside the sphere around the computed center with the radius grown by the
		// center error. The relative margin covers the roundings of these few operations and of the
		// later additions of the origin and the radius. One quantum covers the rounding of the origin.
		r_radius = (Math::sqrt(center_length2) + 2.0 * center_error) * (1.0 + 1e-6) + 1.0;
		for (int k = 0; k < 3; k++) {
			r_center[k] += double(origin[k]);
		}
		return true;
	}

	// Maps a quantized coordinate to the acceleration grid. Points and circumsphere bounds use this same
	// function, it is monotonic, so a point inside the bounds always lands in a cell inside the bounds.
	_FORCE_INLINE_ static int32_t grid_index(double p_coordinate, double p_scale) {
		return int32_t(CLAMP(p_coordinate * p_scale, 0.0, double(ACCEL_GRID_SIZE - 1)));
	}

	static bool simplex_is_coplanar(const Vector3 *p_points, const Simplex &p_simplex) {
		// Checking every possible distance like this is overkill, but only checking
		// one is not enough. If the simplex is almost planar then the vectors p1-p2
		// and p1-p3 can be practically collinear, which makes Plane unreliable.
		for (uint32_t i = 0; i < 4; i++) {
			Plane p(p_points[p_simplex.points[i]], p_points[p_simplex.points[(i + 1) % 4]], p_points[p_simplex.points[(i + 2) % 4]]);
			// This tolerance should not be smaller than the one used with
			// Plane::has_point() when creating the LightmapGI probe BSP tree.
			if (Math::abs(p.distance_to(p_points[p_simplex.points[(i + 3) % 4]])) < 0.001) {
				return true;
			}
		}

		return false;
	}

	static void register_in_grid(List<Simplex *> (&p_grid)[ACCEL_GRID_SIZE][ACCEL_GRID_SIZE][ACCEL_GRID_SIZE], const Vector3i *p_points, const double *p_grid_scale, Simplex *p_simplex) {
		Vector3i from = Vector3i(0, 0, 0);
		Vector3i to = Vector3i(ACCEL_GRID_SIZE - 1, ACCEL_GRID_SIZE - 1, ACCEL_GRID_SIZE - 1);
		double center[3];
		double radius;
		if (circum_sphere_bound(p_points, *p_simplex, center, radius)) {
			for (int k = 0; k < 3; k++) {
				from[k] = grid_index(center[k] - radius, p_grid_scale[k]);
				to[k] = grid_index(center[k] + radius, p_grid_scale[k]);
			}
		}

		for (int32_t x = from.x; x <= to.x; x++) {
			for (int32_t y = from.y; y <= to.y; y++) {
				for (int32_t z = from.z; z <= to.z; z++) {
					GridPos gp;
					gp.pos = Vector3i(x, y, z);
					gp.E = p_grid[x][y][z].push_back(p_simplex);
					p_simplex->grid_positions.push_back(gp);
				}
			}
		}
	}

public:
	struct OutputSimplex {
		uint32_t points[4];
	};

	static Vector<OutputSimplex> tetrahedralize(const Vector<Vector3> &p_points) {
		uint32_t point_count = p_points.size();
		LocalVector<Vector3i> quantized;
		quantized.resize(point_count + 4);
		const Vector3 *src_points = p_points.ptr();
		double grid_scale[3]; // From quantized coordinates to acceleration grid cells.

		{
			AABB rect;
			for (uint32_t i = 0; i < point_count; i++) {
				Vector3 point = src_points[i];
				if (i == 0) {
					rect.position = point;
				} else {
					rect.expand_to(point);
				}
			}

			real_t longest_axis = rect.size[rect.get_longest_axis_index()];
			for (int k = 0; k < 3; k++) {
				// The grid spans the bounds of the points on each axis. A flat axis gets a single cell.
				grid_scale[k] = rect.size[k] > 0 ? double(longest_axis) / double(rect.size[k]) * double(ACCEL_GRID_SIZE) / double(QUANTIZATION_MAX) : 0.0;
			}

			for (uint32_t i = 0; i < point_count; i++) {
				// Scale points to the unit cube and quantize them to integers, so that the
				// predicates can be evaluated exactly.
				quantized[i] = Vector3i((src_points[i] - rect.position) / longest_axis * QUANTIZATION_MAX);
			}

			const int32_t delta_max = int32_t(Math::sqrt(2.0) * 100.0 * double(QUANTIZATION_MAX));
			const int32_t center = QUANTIZATION_MAX / 2;

			// The larger the root simplex is, the more likely it is that the
			// triangulation is convex. If it's not absolutely huge, there can
			// be missing simplices that are not created for the outermost faces
			// of the point cloud if the point density is very low there.
			quantized[point_count + 0] = Vector3i(center, center + delta_max, center);
			quantized[point_count + 1] = Vector3i(center, center - delta_max, center + delta_max);
			quantized[point_count + 2] = Vector3i(center + delta_max, center - delta_max, center - delta_max);
			quantized[point_count + 3] = Vector3i(center - delta_max, center - delta_max, center - delta_max);
		}

		List<Simplex *> acceleration_grid[ACCEL_GRID_SIZE][ACCEL_GRID_SIZE][ACCEL_GRID_SIZE];

		List<Simplex *> simplex_list;
		{
			//create root simplex
			Simplex *root = memnew(Simplex(point_count + 0, point_count + 1, point_count + 2, point_count + 3));
			root->orientation = orientation_sign(quantized.ptr(), *root);
			root->SE = simplex_list.push_back(root);

			for (uint32_t i = 0; i < ACCEL_GRID_SIZE; i++) {
				for (uint32_t j = 0; j < ACCEL_GRID_SIZE; j++) {
					for (uint32_t k = 0; k < ACCEL_GRID_SIZE; k++) {
						GridPos gp;
						gp.E = acceleration_grid[i][j][k].push_back(root);
						gp.pos = Vector3i(i, j, k);
						root->grid_positions.push_back(gp);
					}
				}
			}
		}

		AHashMap<Triangle, uint32_t, TriangleHasher> triangles_inserted;
		LocalVector<Triangle> triangles;

		for (uint32_t i = 0; i < point_count; i++) {
			bool unique = true;
			for (uint32_t j = i + 1; j < point_count; j++) {
				if (quantized[i] == quantized[j]) {
					unique = false;
					break;
				}
			}
			if (!unique) {
				continue;
			}

			Vector3i grid_pos;
			for (int k = 0; k < 3; k++) {
				grid_pos[k] = grid_index(quantized[i][k], grid_scale[k]);
			}

			for (List<Simplex *>::Element *E = acceleration_grid[grid_pos.x][grid_pos.y][grid_pos.z].front(); E;) {
				List<Simplex *>::Element *N = E->next(); //may be deleted

				Simplex *simplex = E->get();

				// Points exactly on the circumsphere are not inside.
				if (insphere_sign(quantized.ptr(), *simplex, i) * simplex->orientation > 0) {
					static const uint32_t triangle_order[4][3] = {
						{ 0, 1, 2 },
						{ 0, 1, 3 },
						{ 0, 2, 3 },
						{ 1, 2, 3 },
					};

					for (uint32_t k = 0; k < 4; k++) {
						Triangle t = Triangle(simplex->points[triangle_order[k][0]], simplex->points[triangle_order[k][1]], simplex->points[triangle_order[k][2]]);
						uint32_t *p = triangles_inserted.getptr(t);
						if (p) {
							// This Delaunay implementation uses the Bowyer-Watson algorithm.
							// The rule is that you don't reuse any triangles that were
							// shared by any of the retriangulated simplices.
							triangles[*p].bad = true;
						} else {
							triangles_inserted.insert(t, triangles.size());
							triangles.push_back(t);
						}
					}

					simplex_list.erase(simplex->SE);

					for (const GridPos &gp : simplex->grid_positions) {
						Vector3i p = gp.pos;
						acceleration_grid[p.x][p.y][p.z].erase(gp.E);
					}
					memdelete(simplex);
				}
				E = N;
			}

			for (const Triangle &triangle : triangles) {
				if (triangle.bad) {
					continue;
				}
				Simplex *new_simplex = memnew(Simplex(triangle.triangle[0], triangle.triangle[1], triangle.triangle[2], i));
				new_simplex->orientation = orientation_sign(quantized.ptr(), *new_simplex);
				// Never flat: the new point can't lie in the plane of a cavity boundary face. The circumspheres
				// of the two simplices sharing the face cut that plane in the same circle, so a point in the
				// plane is either inside both or outside both, and only one of them is in the cavity. This
				// relies on the exact predicates and on the conservative grid, which finds the whole cavity.
				DEV_ASSERT(new_simplex->orientation != 0);
				new_simplex->SE = simplex_list.push_back(new_simplex);
				register_in_grid(acceleration_grid, quantized.ptr(), grid_scale, new_simplex);
			}

			triangles.clear();
			triangles_inserted.clear();
		}

		//print_line("end with simplices: " + itos(simplex_list.size()));
		Vector<OutputSimplex> ret_simplices;
		ret_simplices.resize(simplex_list.size());
		OutputSimplex *ret_simplicesw = ret_simplices.ptrw();
		uint32_t simplices_written = 0;

		for (Simplex *simplex : simplex_list) {
			bool invalid = false;
			for (int j = 0; j < 4; j++) {
				if (simplex->points[j] >= point_count) {
					invalid = true;
					break;
				}
			}
			if (invalid || simplex_is_coplanar(src_points, *simplex)) {
				memdelete(simplex);
				continue;
			}

			ret_simplicesw[simplices_written].points[0] = simplex->points[0];
			ret_simplicesw[simplices_written].points[1] = simplex->points[1];
			ret_simplicesw[simplices_written].points[2] = simplex->points[2];
			ret_simplicesw[simplices_written].points[3] = simplex->points[3];
			simplices_written++;
			memdelete(simplex);
		}

		ret_simplices.resize(simplices_written);

		return ret_simplices;
	}
};
