/*
* LICENCE
* copyright 2014 ~ ****
* Some rights reserved.
* Author: HUFANGYUAN
* Released under CC BY-NC
*/
//Particle_x.h
///defination of class Particle_x

#pragma once
#include "Header.h"
#include "Particle.h"
#include "Polynomial.h"
#include "Derivative.h"
#include <algorithm>
#include <queue>
#include <random>
#define BOOST_PYTHON_STATIC_LIB
#define BOOST_NUMPY_STATIC_LIB
//#define BOOST_ALL_DYN_LINK
//#define BOOST_PYTHON_DYNAMIC_LIB
//#define BOOST_NUMPY_DYNAMIC_LIB
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <Python.h>

	template <typename R, int D, int P>
	class Particle_x : public Particle<R,D> {};

	template <typename R, int P>
	class Particle_x<R,1,P> :  public Particle<R,1>{};

	template <typename R, int P>
	class Particle_x<R,2,P> : public Particle<R,2> {
		typedef mMath::Polynomial_A<R,2,P> PN;
		typedef mMath::Derivative_A<R,2,P> DR;
		typedef mMath::Polynomial_A<R,2,P+1> PNH;
		typedef mMath::Derivative_A<R,2,P+1> DRH;
		typedef Eigen::Matrix<int,2,1> iVec;
		typedef Eigen::Matrix<R,2,1> Vec;
		typedef Eigen::Matrix<R,2,2> Mat;
		typedef Eigen::Matrix<R,PN::value,1> VecP;
		typedef Eigen::Matrix<R,PN::value,2> MatPD;
		typedef Eigen::Matrix<R,PN::value,PN::value> MatPP;
	public:
		Particle_x() : Particle(), python_initialized(false), numpy_initialized(false) {}
		~Particle_x() {}

		__forceinline const R ww(const R& r) const {
				return pow(std::max(1 - r / r0, 0.0), 2);
		}

		__forceinline void poly(const R* in, R* out) const { PN::Run(varrho, in, out); }
		__forceinline void polyH(const R* in, R* out) const { PNH::Run(varrho, in, out); }

		const R DerX(const R* const phi, const int& p) const {
			VecP vv = VecP::Zero();
			const int cx = cell->pos2cell(pos[0][p]);
			const int cy = cell->pos2cell(pos[1][p]);
			for (int i = 0; i < cell->blockSize::value; i++) {
				const int key = cell->hash(cx, cy, i);
				for (int m = 0; m < cell->linkList[key].size(); m++) {
					const int q = cell->linkList[key][m];
					if (type[q] != BD2) {
						const R dr[2] = { pos[0][q] - pos[0][p], pos[1][q] - pos[1][p] };
						const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
						if (dr1 > r0) continue;
						const R w = ww(dr1);
						VecP npq;
						poly(dr, npq.data());
						vv += w * (phi[q] - phi[p]) * npq;
					}
				}
			}
			const VecP aa = invMat[p] * vv;
			return (pn_px_o * aa);
		}

		const R DerY(const R* const phi, const int& p) const {
			VecP vv = VecP::Zero();
			const int cx = cell->pos2cell(pos[0][p]);
			const int cy = cell->pos2cell(pos[1][p]);
			for (int i = 0; i < cell->blockSize::value; i++) {
				const int key = cell->hash(cx, cy, i);
				for (int m = 0; m < cell->linkList[key].size(); m++) {
					const int q = cell->linkList[key][m];
					if (type[q] != BD2) {
						const R dr[2] = { pos[0][q] - pos[0][p], pos[1][q] - pos[1][p] };
						const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
						if (dr1 > r0) continue;
						const R w = ww(dr1);
						VecP npq;
						poly(dr, npq.data());
						vv += w * (phi[q] - phi[p]) * npq;
					}
				}
			}
			const VecP aa = invMat[p] * vv;
			return (pn_py_o * aa);
		}

		const R DerXX(const R* const phi, const int& p) const {
			VecP vv = VecP::Zero();
			const int cx = cell->pos2cell(pos[0][p]);
			const int cy = cell->pos2cell(pos[1][p]);
			for (int i = 0; i < cell->blockSize::value; i++) {
				const int key = cell->hash(cx, cy, i);
				for (int m = 0; m < cell->linkList[key].size(); m++) {
					const int q = cell->linkList[key][m];
					if (type[q] != BD2) {
						const R dr[2] = { pos[0][q] - pos[0][p], pos[1][q] - pos[1][p] };
						const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
						if (dr1 > r0) continue;
						const R w = ww(dr1);
						VecP npq;
						poly(dr, npq.data());
						vv += w * (phi[q] - phi[p]) * npq;
					}
				}
			}
			const VecP aa = invMat[p] * vv;
			return (pn_pxx_o * aa);
		}

		const R DerYY(const R* const phi, const int& p) const {
			VecP vv = VecP::Zero();
			const int cx = cell->pos2cell(pos[0][p]);
			const int cy = cell->pos2cell(pos[1][p]);
			for (int i = 0; i < cell->blockSize::value; i++) {
				const int key = cell->hash(cx, cy, i);
				for (int m = 0; m < cell->linkList[key].size(); m++) {
					const int q = cell->linkList[key][m];
					if (type[q] != BD2) {
						const R dr[2] = { pos[0][q] - pos[0][p], pos[1][q] - pos[1][p] };
						const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
						if (dr1 > r0) continue;
						const R w = ww(dr1);
						VecP npq;
						poly(dr, npq.data());
						vv += w * (phi[q] - phi[p]) * npq;
					}
				}
			}
			const VecP aa = invMat[p] * vv;
			return (pn_pyy_o * aa);
		}

		const Vec Grad(const R* const phi, const int& p) const {
			VecP vv = VecP::Zero();
			const int cx = cell->pos2cell(pos[0][p]);
			const int cy = cell->pos2cell(pos[1][p]);
			for (int i = 0; i < cell->blockSize::value; i++) {
				const int key = cell->hash(cx, cy, i);
				for (int m = 0; m < cell->linkList[key].size(); m++) {
					const int q = cell->linkList[key][m];
					if (type[q] != BD2) {
						const R dr[2] = { pos[0][q] - pos[0][p], pos[1][q] - pos[1][p] };
						const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
						if (dr1 > r0) continue;
						const R w = ww(dr1);
						VecP npq;
						poly(dr, npq.data());
						vv += w * (phi[q] - phi[p]) * npq;
					}
				}
			}
			const VecP aa = invMat[p] * vv;
			return (pn_p_o * aa);
		}

		const Vec Grad(const R* const phi, const int& p, const int& mask) const {
			MatPP mm = MatPP::Zero();
			VecP vv = VecP::Zero();
			const int cx = cell->pos2cell(pos[0][p]);
			const int cy = cell->pos2cell(pos[1][p]);
			for (int i = 0; i < cell->blockSize::value; i++) {
				const int key = cell->hash(cx, cy, i);
				for (int m = 0; m < cell->linkList[key].size(); m++) {
					const int q = cell->linkList[key][m];
					if (type[q] & mask) {
						const R dr[2] = { pos[0][q] - pos[0][p], pos[1][q] - pos[1][p] };
						const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
						if (dr1 > r0) continue;
						const R w = ww(dr1);
						VecP npq;
						poly(dr, npq.data());
						vv += (w * (phi[q] - phi[p])) * npq;
						mm += (w* npq) * npq.transpose();
					}
				}
			}
			MatPP inv = MatPP::Zero();
			if (abs(mm.determinant()) < eps_mat) {
#if DEBUG
				std::cout << " Determinant defficiency: "; PRINT(p);
#endif
				auto mm_ = mm.block<2, 2>(0, 0);
				if (abs(mm_.determinant()) < eps_mat) inv = MatPP::Zero();
				else inv.block<2, 2>(0, 0) = mm_.inverse();
			}
			else inv = mm.inverse();
			const VecP aa = inv * vv;
			return (pn_p_o * aa);
		}

		const R Div(const R* const phix, const R* const phiy, const int& p) const {
			VecP vvx = VecP::Zero();
			VecP vvy = VecP::Zero();
			const int cx = cell->pos2cell(pos[0][p]);
			const int cy = cell->pos2cell(pos[1][p]);
			for (int i = 0; i < cell->blockSize::value; i++) {
				const int key = cell->hash(cx, cy, i);
				for (int m = 0; m < cell->linkList[key].size(); m++) {
					const int q = cell->linkList[key][m];
					if (type[q] != BD2) {
						const R dr[2] = { pos[0][q] - pos[0][p], pos[1][q] - pos[1][p] };
						const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
						if (dr1 > r0) continue;
						const R w = ww(dr1);
						VecP npq;
						poly(dr, npq.data());
						vvx += w * (phix[q] - phix[p]) * npq;
						vvy += w * (phiy[q] - phiy[p]) * npq;
					}
				}
			}
			const R pupx = pn_px_o* invMat[p] * vvx;
			const R pvpy = pn_py_o* invMat[p] * vvy;
			return pupx + pvpy;
		}

		const R Lap(const R* const phi, const int& p) const {
			VecP vv = VecP::Zero();
			const int cx = cell->pos2cell(pos[0][p]);
			const int cy = cell->pos2cell(pos[1][p]);
			for (int i = 0; i < cell->blockSize::value; i++) {
				const int key = cell->hash(cx, cy, i);
				for (int m = 0; m < cell->linkList[key].size(); m++) {
					const int q = cell->linkList[key][m];
					if (type[q] != BD2) {
						const R dr[2] = { pos[0][q] - pos[0][p], pos[1][q] - pos[1][p] };
						const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
						if (dr1 > r0) continue;
						const R w = ww(dr1);
						VecP npq;
						poly(dr, npq.data());
						vv += w * (phi[q] - phi[p]) * npq;
					}
				}
			}
			const VecP aa = invMat[p] * vv;
			return (pn_lap_o * aa);
		}

		const Vec Lap(const R* const phix, const R* const phiy, const int& p) const {
			VecP vvx = VecP::Zero();
			VecP vvy = VecP::Zero();
			const int cx = cell->pos2cell(pos[0][p]);
			const int cy = cell->pos2cell(pos[1][p]);
			for (int i = 0; i < cell->blockSize::value; i++) {
				const int key = cell->hash(cx, cy, i);
				for (int m = 0; m < cell->linkList[key].size(); m++) {
					const int q = cell->linkList[key][m];
					if (type[q] != BD2) {
						const R dr[2] = { pos[0][q] - pos[0][p], pos[1][q] - pos[1][p] };
						const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
						if (dr1 > r0) continue;
						const R w = ww(dr1);
						VecP npq;
						poly(dr, npq.data());
						vvx += w * (phix[q] - phix[p]) * npq;
						vvy += w * (phiy[q] - phiy[p]) * npq;
					}
				}
			}
			const VecP aax = invMat[p] * vvx;
			const VecP aay = invMat[p] * vvy;
			Vec ret;
			ret[0] = pn_lap_o * aax;
			ret[1] = pn_lap_o * aay;
			return ret;
		}

		const R Rot(const R* const phix, const R* const phiy, const int& p) const {
			VecP vvx = VecP::Zero();
			VecP vvy = VecP::Zero();
			const int cx = cell->pos2cell(pos[0][p]);
			const int cy = cell->pos2cell(pos[1][p]);
			for (int i = 0; i < cell->blockSize::value; i++) {
				const int key = cell->hash(cx, cy, i);
				for (int m = 0; m < cell->linkList[key].size(); m++) {
					const int q = cell->linkList[key][m];
					if (type[q] != BD2) {
						const R dr[2] = { pos[0][q] - pos[0][p], pos[1][q] - pos[1][p] };
						const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
						if (dr1 > r0) continue;
						const R w = ww(dr1);
						VecP npq;
						poly(dr, npq.data());
						vvx += w * (phix[q] - phix[p]) * npq;
						vvy += w * (phiy[q] - phiy[p]) * npq;
					}
				}
			}
			const R pupy = pn_py_o* invMat[p] * vvx;
			const R pvpx = pn_px_o* invMat[p] * vvy;
			return (pupy - pvpx);
		}

		const R interpolateLSA(const R* const phi, const R& px, const R& py) const {
			int p = -1;
			R dp2_min = 2 * 2 * dp*dp;
			R dis2 = std::numeric_limits<R>::max();
			const int cx = cell->pos2cell(px);
			const int cy = cell->pos2cell(py);
			for (int i = 0; i < cell->blockSize::value; i++) {
				const int key = cell->hash(cx, cy, i);
				for (int m = 0; m < cell->linkList[key].size(); m++) {
					const int q = cell->linkList[key][m];
					if (type[q] != BD2) {
						const R dr[2] = { pos[0][q] - px, pos[1][q] - py };
						const R dr2 = (dr[0] * dr[0] + dr[1] * dr[1]);
						const R dr1 = sqrt(dr2);
						if (dr1 > r0) continue;
						if (dr2 < dis2) {
							dis2 = dr2;
							p = q;
						}
					}
				}
			}
			if (p == -1 || type[p] == BD2 || dis2 > dp2_min) return R(0);
			return interpolateLSA(phi, p, px, py);
		}

		const Vec interpolateLSA(const R* const phix, const R* const phiy, const R& px, const R& py) const {
			int p = -1;
			R dp2_min = 2 * 2 * dp*dp;
			R dis2 = std::numeric_limits<R>::max();
			const int cx = cell->pos2cell(px);
			const int cy = cell->pos2cell(py);
			for (int i = 0; i < cell->blockSize::value; i++) {
				const int key = cell->hash(cx, cy, i);
				for (int m = 0; m < cell->linkList[key].size(); m++) {
					const int q = cell->linkList[key][m];
					const R dr[2] = { pos[0][q] - px, pos[1][q] - py };
					const R dr2 = (dr[0] * dr[0] + dr[1] * dr[1]);
					const R dr1 = sqrt(dr2);
					if (dr1 > r0) continue;
					if (dr2 < dis2) {
						dis2 = dr2;
						p = q;
					}
				}
			}
			if (p == -1 || type[p]== BD2 || dis2 > dp2_min) return Vec::Zero();
			return interpolateLSA(phix, phiy, p, px, py);
		}

		const R interpolateLSA(const R* const phi, const int& p, const R& px, const R& py) const {
			const R dx = px - pos[0][p];
			const R dy = py - pos[1][p];
			MatPP mm = MatPP::Zero();
			VecP vv = VecP::Zero();
			const int cx = cell->pos2cell(pos[0][p]);
			const int cy = cell->pos2cell(pos[1][p]);
			for (int i = 0; i < cell->blockSize::value; i++) {
				const int key = cell->hash(cx, cy, i);
				for (int m = 0; m < cell->linkList[key].size(); m++) {
					const int q = cell->linkList[key][m];
					if (type[q] != BD2) {
						const R dr[2] = { pos[0][q] - pos[0][p], pos[1][q] - pos[1][p] };
						const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
						if (dr1 > r0) continue;
						const R w = ww(dr1);
						VecP npq;
						poly(dr, npq.data());
						mm += (w * npq) * npq.transpose();
						vv += w * (phi[q] - phi[p]) * npq;
					}
				}
			}
			MatPP inv = MatPP::Zero();
			if (abs(mm.determinant()) < eps_mat) {
#if DEBUG
				std::cout << " (interpolateLSA) ID: " << p << " --- " << " Determinant defficiency: " << mm.determinant() << std::endl;
#endif
				auto mm_ = mm.block<2, 2>(0, 0);
				if (abs(mm_.determinant()) < eps_mat) {
					return phi[p];
				}
				else inv.block<2, 2>(0, 0) = mm_.inverse();
			}
			else inv = mm.inverse();

			const VecP aa = inv * vv;
			const R Px = pn_px_o* aa;
			const R Py = pn_py_o* aa;
			const R Pxx = pn_pxx_o* aa;
			const R Pxy = pn_pxy_o* aa;
			const R Pyy = pn_pyy_o* aa;
			return phi[p] + (dx*Px + dy*Py) + 0.5* (dx*dx*Pxx + R(2)*dx*dy*Pxy + dy*dy*Pyy);
		}

		const Vec interpolateLSA(const R* const phix, const R* const phiy, const int& p, const R& px, const R& py) const {
			const R dx = px - pos[0][p];
			const R dy = py - pos[1][p];
			MatPP mm = MatPP::Zero();
			VecP vvx = VecP::Zero();
			VecP vvy = VecP::Zero();
			const int cx = cell->pos2cell(pos[0][p]);
			const int cy = cell->pos2cell(pos[1][p]);
			for (int i = 0; i < cell->blockSize::value; i++) {
				const int key = cell->hash(cx, cy, i);
				for (int m = 0; m < cell->linkList[key].size(); m++) {
					const int q = cell->linkList[key][m];
					if (type[q] != BD2) {
						const R dr[2] = { pos[0][q] - pos[0][p], pos[1][q] - pos[1][p] };
						const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
						if (dr1 > r0) continue;
						const R w = ww(dr1);
						VecP npq;
						poly(dr, npq.data());
						mm += (w * npq) * npq.transpose();
						vvx += w * (phix[q] - phix[p]) * npq;
						vvy += w * (phiy[q] - phiy[p]) * npq;
					}
				}
			}
			MatPP inv = MatPP::Zero();
			if (abs(mm.determinant()) < eps_mat) {
#if DEBUG
				std::cout << " (interpolateLSA) ID: " << p << " --- " << " Determinant defficiency: " << mm.determinant() << std::endl;
#endif
				auto mm_ = mm.block<2, 2>(0, 0);
				if (abs(mm_.determinant()) < eps_mat) {
					Vec retv;
					retv[0] = phix[p];
					retv[1] = phiy[p];
					return retv;
				}
				else inv.block<2, 2>(0, 0) = mm_.inverse();
			}
			else inv = mm.inverse();

			const VecP aax = inv * vvx;
			const VecP aay = inv * vvy;
			const R Px[2] = { pn_px_o* aax, pn_px_o* aay };
			const R Py[2] = { pn_py_o* aax, pn_py_o* aay };
			const R Pxx[2] = { pn_pxx_o* aax, pn_pxx_o* aay };
			const R Pxy[2] = { pn_pxy_o* aax, pn_pxy_o* aay };
			const R Pyy[2] = { pn_pyy_o* aax, pn_pyy_o* aay };
			Vec ret;
			ret[0] = phix[p] + (dx*Px[0] + dy*Py[0]) + 0.5* (dx*dx*Pxx[0] + R(2)*dx*dy*Pxy[0] + dy*dy*Pyy[0]);
			ret[1] = phiy[p] + (dx*Px[1] + dy*Py[1]) + 0.5* (dx*dx*Pxx[1] + R(2)*dx*dy*Pxy[1] + dy*dy*Pyy[1]);
			return ret;
		}

		const R interpolateLSAU(const R* const phi, const int& p, const R& px, const R& py) const {
			const R dx = px - pos[0][p];
			const R dy = py - pos[1][p];
			MatPP mm = MatPP::Zero();
			VecP vv = VecP::Zero();
			const int cx = cell->pos2cell(pos[0][p]);
			const int cy = cell->pos2cell(pos[1][p]);
			for (int i = 0; i < cell->blockSize::value; i++) {
				const int key = cell->hash(cx, cy, i);
				for (int m = 0; m < cell->linkList[key].size(); m++) {
					const int q = cell->linkList[key][m];
					if (type[q] != BD2) {
						const R dr[2] = { pos[0][q] - pos[0][p], pos[1][q] - pos[1][p] };
						const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
						if (dr1 > r0 || (dx*dr[0] + dy*dr[1]) < R(0)) continue;
						const R w = ww(dr1);
						VecP npq;
						poly(dr, npq.data());
						mm += (w * npq) * npq.transpose();
						vv += w * (phi[q] - phi[p]) * npq;
					}
				}
			}
			MatPP inv = MatPP::Zero();
			if (abs(mm.determinant()) < eps_mat) {
#if DEBUG
				std::cout << " (interpolateLSAU) ID: " << p << " --- " << " Determinant defficiency: " << mm.determinant() << std::endl;
#endif
				auto mm_ = mm.block<2, 2>(0, 0);
				if (abs(mm_.determinant()) < eps_mat) {
					return phi[p];
				}
				else inv.block<2, 2>(0, 0) = mm_.inverse();
			}
			else inv = mm.inverse();

			const VecP aa = inv * vv;
			const R Px = pn_px_o* aa;
			const R Py = pn_py_o* aa;
			const R Pxx = pn_pxx_o* aa;
			const R Pxy = pn_pxy_o* aa;
			const R Pyy = pn_pyy_o* aa;
			return phi[p] + (dx*Px + dy*Py) + 0.5* (dx*dx*Pxx + R(2)*dx*dy*Pxy + dy*dy*Pyy);
		}

		const Vec interpolateLSAU(const R* const phix, const R* const phiy, const int& p, const R& px, const R& py) const {
			const R dx = px - pos[0][p];
			const R dy = py - pos[1][p];
			MatPP mm = MatPP::Zero();
			VecP vvx = VecP::Zero();
			VecP vvy = VecP::Zero();
			const int cx = cell->pos2cell(pos[0][p]);
			const int cy = cell->pos2cell(pos[1][p]);
			for (int i = 0; i < cell->blockSize::value; i++) {
				const int key = cell->hash(cx, cy, i);
				for (int m = 0; m < cell->linkList[key].size(); m++) {
					const int q = cell->linkList[key][m];
					if (type[q] != BD2) {
						const R dr[2] = { pos[0][q] - pos[0][p], pos[1][q] - pos[1][p] };
						const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
						if (dr1 > r0 || (dx*dr[0] + dy*dr[1]) < R(0)) continue;
						const R w = ww(dr1);
						VecP npq;
						poly(dr, npq.data());
						mm += (w * npq) * npq.transpose();
						vvx += w * (phix[q] - phix[p]) * npq;
						vvy += w * (phiy[q] - phiy[p]) * npq;
					}
				}
			}
			MatPP inv = MatPP::Zero();
			if (abs(mm.determinant()) < eps_mat) {
#if DEBUG
				std::cout << " (interpolateLSAU) ID: " << p << " --- " << " Determinant defficiency: " << mm.determinant() << std::endl;
#endif
				auto mm_ = mm.block<2, 2>(0, 0);
				if (abs(mm_.determinant()) < eps_mat) {
					Vec retv;
					retv[0] = phix[p];
					retv[1] = phiy[p];
					return retv;
				}
				else inv.block<2, 2>(0, 0) = mm_.inverse();
			}
			else inv = mm.inverse();

			const VecP aax = inv * vvx;
			const VecP aay = inv * vvy;
			const R Px[2] = { pn_px_o* aax, pn_px_o* aay };
			const R Py[2] = { pn_py_o* aax, pn_py_o* aay };
			const R Pxx[2] = { pn_pxx_o* aax, pn_pxx_o* aay };
			const R Pxy[2] = { pn_pxy_o* aax, pn_pxy_o* aay };
			const R Pyy[2] = { pn_pyy_o* aax, pn_pyy_o* aay };
			Vec ret;
			ret[0] = phix[p] + (dx*Px[0] + dy*Py[0]) + 0.5* (dx*dx*Pxx[0] + R(2)*dx*dy*Pxy[0] + dy*dy*Pyy[0]);
			ret[1] = phiy[p] + (dx*Px[1] + dy*Py[1]) + 0.5* (dx*dx*Pxx[1] + R(2)*dx*dy*Pxy[1] + dy*dy*Pyy[1]);
			return ret;
		}

		template <int StencilsX = 1, int StencilsY = 3, int Stencils = StencilsX*StencilsY>
		const R interpolateWENO_A_(const R* const phi, const int& p, const R& px, const R& py) const {
			const R dx = px - pos[0][p];
			const R dy = py - pos[1][p];
			const R dd = sqrt(dx*dx + dy*dy);
			if (dd < eps) return phi[p];
			const R upx = dx / dd;
			const R upy = dy / dd;
			const R alpha = 2.* M_PI / StencilsX;
			Vec dir[StencilsX];
			Vec ctr[Stencils];
			for (int i = 0; i < StencilsX; i++) {
				const R theta = i* alpha;
				const R ct = cos(theta);
				const R st = sin(theta);
				dir[i] << ct*upx + st*upy, ct*upx - st*upy;
			}
			for (int j = 0; j < StencilsY; j++) {
				//const R dis = r0* ( R(1.) - R(2.)*(j + 1) / (1 + StencilsY) );
				const R dis = r0* (R(1.) - R(1.)*(j + 1) / (StencilsY));
				for (int i = 0; i < StencilsX; i++) {
					const int stcId = i* StencilsY + j;
					ctr[stcId][0] = pos[0][p] + dis*dir[i][0];
					ctr[stcId][1] = pos[1][p] + dis*dir[i][1];
				}
			}
			MatPP mm[Stencils];
			VecP vv[Stencils];
			for (int i = 0; i < Stencils; i++) {
				mm[i] = MatPP::Zero();
				vv[i] = VecP::Zero();
			}
			const int cx = cell->pos2cell(pos[0][p]);
			const int cy = cell->pos2cell(pos[1][p]);
			for (int i = 0; i < cell->blockSize::value; i++) {
				const int key = cell->hash(cx, cy, i);
				for (int m = 0; m < cell->linkList[key].size(); m++) {
					const int q = cell->linkList[key][m];
					if (type[q] == BD2) continue;
					const R dr[2] = { pos[0][q] - pos[0][p], pos[1][q] - pos[1][p] };
					const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
					if (dr1 > r0) continue;
					for (int stcId = 0; stcId < Stencils; stcId++) {
						const R disx = (pos[0][q] - ctr[stcId][0]);
						const R disy = (pos[1][q] - ctr[stcId][1]);
						const R dis = sqrt(disx*disx + disy*disy);
						const R w = ww(dis);
						VecP npq;
						part->poly(dr, npq.data());
						mm[stcId] += (w* npq)* npq.transpose();
						vv[stcId] += (w* npq)* (phi[q] - phi[p]);
					}
				}
			}
			R oscillationIndicator[Stencils];
			R stencilWeight[Stencils];
			R stencilWeightNorm[Stencils];
			VecP polyCoef[Stencils];
			for (int i = 0; i < Stencils; i++) {
				MatPP inv = MatPP::Zero();
				if (abs(mm[i].determinant()) < part->eps_mat) {
					auto mm_ = mm[i].block<2, 2>(0, 0);
					if (abs(mm_.determinant()) < part->eps_mat) {
						inv = MatPP::Zero();
						inv(0, 0) = R(1) / eps;
						inv(1, 1) = R(1) / eps;
					}
					else inv.block<2, 2>(0, 0) = mm_.inverse();
				}
				else inv = mm[i].inverse();
				polyCoef[i] = inv * vv[i];
				const int offset = PN::value - mMath::H<D, P>::value;
				oscillationIndicator[i] = R(0.);
				for (int term = offset; term < PN::value; term++) {
					oscillationIndicator[i] += abs(polyCoef[i][term]);
				}
			}
			const R epsilon = 1.e-6;
			const int magnifier = 5;
			for (int i = 0; i < Stencils; i++) {
				stencilWeight[i] = R(1) / pow(epsilon + oscillationIndicator[i], magnifier);
			}
			R stencilWeightSum = R(0);
			for (int i = 0; i < Stencils; i++) {
				stencilWeightSum += stencilWeight[i];
			}
			for (int i = 0; i < Stencils; i++) {
				stencilWeightNorm[i] = stencilWeight[i] / stencilWeightSum;
			}
			VecP combinedCoef = VecP::Zero();
			for (int i = 0; i < Stencils; i++) {
				combinedCoef += stencilWeightNorm[i] * polyCoef[i];
			}
			const R Px = pn_px_o* combinedCoef;
			const R Py = pn_py_o* combinedCoef;
			const R Pxx = pn_pxx_o* combinedCoef;
			const R Pxy = pn_pxy_o* combinedCoef;
			const R Pyy = pn_pyy_o* combinedCoef;
			return phi[p] + (dx*Px + dy*Py) + 0.5* (dx*dx*Pxx + 2.0*dx*dy*Pxy + dy*dy*Pyy);
		}

		template <int StencilsX = 1, int StencilsY = 3, int Stencils = StencilsX*StencilsY>
		const R interpolateWENO_B_(const R* const phi, const int& p, const R& px, const R& py) const {
			const R dx = px - pos[0][p];
			const R dy = py - pos[1][p];
			const R dd = sqrt(dx*dx + dy*dy);
			if (dd < eps) return phi[p];
			const R upx = dx / dd;
			const R upy = dy / dd;
			const R alpha = 2.* M_PI / StencilsX;
			Vec dir[StencilsX];
			Vec ctr[Stencils];
			for (int i = 0; i < StencilsX; i++) {
				const R theta = i* alpha;
				const R ct = cos(theta);
				const R st = sin(theta);
				dir[i] << ct*upx + st*upy, ct*upx - st*upy;
			}
			for (int j = 0; j < StencilsY; j++) {
				//const R dis = r0* ( R(1.) - R(2.)*(j + 1) / (1 + StencilsY) );
				const R dis = r0* (R(1.) - R(1.)*(j + 1) / (StencilsY));
				for (int i = 0; i < StencilsX; i++) {
					const int stcId = i* StencilsY + j;
					ctr[stcId][0] = pos[0][p] + dis*dir[i][0];
					ctr[stcId][1] = pos[1][p] + dis*dir[i][1];
				}
			}
			MatPP mm[Stencils];
			VecP vv[Stencils];
			for (int i = 0; i < Stencils; i++) {
				mm[i] = MatPP::Zero();
				vv[i] = VecP::Zero();
			}
			const int cx = cell->pos2cell(pos[0][p]);
			const int cy = cell->pos2cell(pos[1][p]);
			for (int i = 0; i < cell->blockSize::value; i++) {
				const int key = cell->hash(cx, cy, i);
				for (int m = 0; m < cell->linkList[key].size(); m++) {
					const int q = cell->linkList[key][m];
					if (type[q] == BD2) continue;
					const R dr[2] = { pos[0][q] - pos[0][p], pos[1][q] - pos[1][p] };
					const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
					if (dr1 > r0) continue;
					for (int stcId = 0; stcId < Stencils; stcId++) {
						const R disx = (pos[0][q] - ctr[stcId][0]);
						const R disy = (pos[1][q] - ctr[stcId][1]);
						const R dis = sqrt(disx*disx + disy*disy);
						const R w = ww(dis);
						VecP npq;
						poly(dr, npq.data());
						mm[stcId] += (w* npq)* npq.transpose();
						vv[stcId] += (w* npq)* (phi[q] - phi[p]);
					}
				}
			}
			R oscillationIndicator[Stencils];
			R stencilWeight[Stencils];
			R stencilWeightNorm[Stencils];
			VecP polyCoef[Stencils];
			for (int i = 0; i < Stencils; i++) {
				MatPP inv = MatPP::Zero();
				if (abs(mm[i].determinant()) < eps_mat) {
					return phi[p];
					auto mm_ = mm[i].block<2, 2>(0, 0);
					if (abs(mm_.determinant()) < eps_mat) {
						inv = MatPP::Zero();
						inv(0, 0) = R(1) / eps;
						inv(1, 1) = R(1) / eps;
					}
					else inv.block<2, 2>(0, 0) = mm_.inverse();
				}
				else inv = mm[i].inverse();
				polyCoef[i] = inv * vv[i];
				oscillationIndicator[i] = R(0.);
			}

			for (auto i = 0; i < Stencils; i++) {
				const R A = polyCoef[i][0] * polyCoef[i][0];
				const R B = polyCoef[i][1] * polyCoef[i][1];
				const R C = polyCoef[i][2] * polyCoef[i][2];
				const R D = polyCoef[i][3] * polyCoef[i][3];
				const R E = polyCoef[i][4] * polyCoef[i][4];
				const R beta1 = (dp*dp*dp)*(A + B) + (dp*dp*dp*dp*dp / R(6))*(R(2)*C + D + R(2)*E);
				const R beta2 = (dp*dp*dp*dp*dp)*(R(4)*C + D + R(4)*E);
				//oscillationIndicator[i] = R(4)*beta1 - (R(1) / R(3))*beta2;
				oscillationIndicator[i] = beta1 + beta2;
			}
			const R epsilon = 1.e-6;
			const int magnifier = 5;
			for (int i = 0; i < Stencils; i++) {
				stencilWeight[i] = R(1) / pow(epsilon + oscillationIndicator[i], magnifier);
			}
			R stencilWeightSum = R(0);
			for (int i = 0; i < Stencils; i++) {
				stencilWeightSum += stencilWeight[i];
			}
			for (int i = 0; i < Stencils; i++) {
				stencilWeightNorm[i] = stencilWeight[i] / stencilWeightSum;
			}
			VecP combinedCoef = VecP::Zero();
			for (int i = 0; i < Stencils; i++) {
				combinedCoef += stencilWeightNorm[i] * polyCoef[i];
			}
			const R Px = pn_px_o* combinedCoef;
			const R Py = pn_py_o* combinedCoef;
			const R Pxx = pn_pxx_o* combinedCoef;
			const R Pxy = pn_pxy_o* combinedCoef;
			const R Pyy = pn_pyy_o* combinedCoef;
			return phi[p] + (dx*Px + dy*Py) + 0.5* (dx*dx*Pxx + 2.0*dx*dy*Pxy + dy*dy*Pyy);
		}

		__forceinline const R interpolateWENO(const R* const phi, const int& p, const R& px, const R& py) const {
			return interpolateWENO_B_<>(phi, p, px, py);
		}

		void updateInvMat() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(np); p++) {
				if (type[p] == BD2) continue;
				MatPP mm = MatPP::Zero();
				const int cx = cell->pos2cell(pos[0][p]);
				const int cy = cell->pos2cell(pos[1][p]);
				for (int i = 0; i < cell->blockSize::value; i++) {
					const int key = cell->hash(cx, cy, i);
					for (int m = 0; m < cell->linkList[key].size(); m++) {
						const int q = cell->linkList[key][m];
						if (type[q] == BD2) continue;
						const R dr[2] = { pos[0][q] - pos[0][p], pos[1][q] - pos[1][p] };
						const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
						if (dr1 > r0) continue;
						const R w = ww(dr1);
						VecP npq;
						poly(dr, npq.data());
						mm += (w* npq) * npq.transpose();
					}
				}
				MatPP& invRef = invMat.at(p);
				invRef = MatPP::Zero();
				if (abs(mm.determinant()) < eps_mat) {
#if DEBUG
					std::cout << " Determinant defficiency: "; PRINT(p);
#endif
					auto mm_ = mm.block<2, 2>(0, 0);
					if (abs(mm_.determinant()) < eps_mat) invRef = MatPP::Zero();
					else invRef.block<2, 2>(0, 0) = mm_.inverse();
				}
				else invRef = mm.inverse();
			}
		}

		int isSurf(int p) {
			if (type[p] == BD2) return 0;
			Mat mm = Mat::Zero();
			//vec gc = vec(0., 0., 0.);
			const int cx = cell->pos2cell(pos[0][p]);
			const int cy = cell->pos2cell(pos[1][p]);
			for (int i = 0; i < cell->blockSize::value; i++) {
				const int key = cell->hash(cx, cy, i);
				for (int m = 0; m < cell->linkList[key].size(); m++) {
					const int q = cell->linkList[key][m];
					if (q == p) continue;
					const R dr[2] = { pos[0][q] - pos[0][p], pos[1][q] - pos[1][p] };
					const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
					if (dr1 > r0) continue;
					const R  w = ww(dr1);
					Vec npq;
					npq[0] = dr[0] / dr1;
					npq[1] = dr[1] / dr1;
					mm += (w* npq)* npq.transpose();
					//gc += w* npq;
				}
			}
			//gc = gc.norm();
			//mm = (R(2) / n0)* mm;
			Eigen::SelfAdjointEigenSolver<Mat> sol(mm);
			Vec eig = sol.eigenvalues();
			Mat eigvs = sol.eigenvectors();
			R eigmin = eig[0];
			Vec neigv1 = eigvs.col(0);
			if (eigmin <= 0.2) return 1;
			if (eigmin > 0.8) return 0;
			Vec neigv2 = -neigv1;
			//if (neigv * gc > 0.) neigv = -1.*neigv;
			//gc = -1. * gc;
			const R root2 = 1.415;
			int flag1 = 1, flag2 = 1;
			for (int i = 0; i < cell->blockSize::value; i++) {
				const int key = cell->hash(cx, cy, i);
				for (int m = 0; m < cell->linkList[key].size(); m++) {
					const int q = cell->linkList[key][m];
					if (q == p) continue;
					const R dr[2] = { pos[0][q] - pos[0][p], pos[1][q] - pos[1][p] };
					const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
					Vec drv, posp, posq;
					drv[0] = dr[0]; drv[1] = dr[1];
					posp[0] = pos[0][p]; posp[1] = pos[1][p];
					posq[0] = pos[0][q]; posq[1] = pos[1][q];
					if (dr1 < root2 * dp) {
						if ((drv / dr1).dot(neigv1) >(root2 / 2.)) flag1 = 0;
						if ((drv / dr1).dot(neigv2) >(root2 / 2.)) flag2 = 0;
					}
					else {
						if ((posp + dp * neigv1 - posq).norm() < dp) flag1 = 0;
						if ((posp + dp * neigv2 - posq).norm() < dp) flag2 = 0;
					}
				}
			}
			if (flag1) {
				//fsNorm[p] = neigv1;
				return 1;
			}
			if (flag2) {
				//fsNorm[p] = neigv2;
				return 1;
			}
			return 0;
		}
		int isSurfML(boost::python::object& global, int p) {
			namespace PY = boost::python;
			namespace NPY = boost::python::numpy;
			static const int N = 8;
			std::vector<int> nbr;
			nNearestNeighbor<N>(nbr, p);
			if (!numpy_initialized) {
				NPY::initialize();
				numpy_initialized = true;
			}
			NPY::ndarray x = NPY::zeros( PY::make_tuple(2*N+1, 1), NPY::dtype::get_builtin<float>() );
			x[0][0] = type[p];
			for (size_t i = 0; i < nbr.size(); i++) {
				x[i * 2 + 1][0] = (pos[0][nbr[i]] - pos[0][p]) / dp;
				x[i * 2 + 2][0] = (pos[1][nbr[i]] - pos[1][p]) / dp;
			}
			for (size_t i = nbr.size(); i < N; i++) {
				x[i * 2 + 1][0] = 0;
				x[i * 2 + 2][0] = 0;
			}
			PY::object nn = global["nn"];
			PY::object predict01 = nn.attr("predict01");
			PY::object ret = predict01(x);
			PY::object np = global["numpy"];
			PY::object asscalar = np.attr("asscalar");
			PY::object iret = asscalar(ret);
			return int( PY::extract<int>(iret) );
		}

		void makeSurf() {
			using namespace boost::python;
			if (!python_initialized) {
				Py_Initialize();
				python_initialized = true;
			}
			object main_module = import("__main__");
			object global = main_module.attr("__dict__");
			exec("import sys", global, global);
			exec("sys.path.append('.\\python')", global, global);
			exec("import numpy", global, global);
			exec("from NN import NN as NN", global, global);
			exec("Layers = (17, 8, 8, 1)", global, global);
			exec("nn = NN(Layers = Layers)", global, global);
			exec("nn.load('.\\python\\config')", global, global);
			for (int p = 0; p < int(np); p++) surf[p] = double(isSurfML(global, p));
		}

		void Redistribute() {
			using namespace boost::python;
			namespace NPY = boost::python::numpy;
			if (!python_initialized) {
				Py_Initialize();
				python_initialized = true;
			}
			if (!numpy_initialized) {
				NPY::initialize();
				numpy_initialized = true;
			}
			object main_module = import("__main__");
			object global = main_module.attr("__dict__");
			exec("import tensorflow as tf", global, global);
			exec("saver = tf.train.Saver()", global, global);
			exec("sess = tf.Session()", global, global);
			exec("saver.restore(sess, '.\\python\\tf_model\\model.ckpt')", global, global);
			const int N = 24;
			std::vector<int> nbr;
			nNearestNeighbor<N>(nbr, p);
			if (!numpy_initialized) {
				NPY::initialize();
				numpy_initialized = true;
			}
			NPY::ndarray x = NPY::zeros(PY::make_tuple(1, 2 * N), NPY::dtype::get_builtin<float>());
			for (size_t i = 0; i < nbr.size(); i++) {
				x[0][i * 2] = (pos[0][nbr[i]] - pos[0][p]) / dp;
				x[0][i * 2 + 1] = (pos[1][nbr[i]] - pos[1][p]) / dp;
			}
			for (size_t i = nbr.size(); i < N; i++) {
				x[0][i * 2] = 0;
				x[0][i * 2 + 1] = 0;
			}
			///predict here
		}

		template <int N = 8>
		void nNearestNeighbor(std::vector<int>& nbr, const int p) const {
			typedef std::pair<R,int> R2i;
			std::priority_queue<R2i, std::vector<R2i>, std::greater<R2i>> que;
			const int cx = cell->pos2cell(pos[0][p]);
			const int cy = cell->pos2cell(pos[1][p]);
			for (int i = 0; i < cell->blockSize::value; i++) {
				const int key = cell->hash(cx, cy, i);
				for (int m = 0; m < cell->linkList[key].size(); m++) {
					const int q = cell->linkList[key][m];
					if (type[q] == BD2) continue;
					const R dr[2] = { pos[0][q] - pos[0][p], pos[1][q] - pos[1][p] };
					const R dr2 = dr[0] * dr[0] + dr[1] * dr[1];
					que.push({ dr2, q });
				}
			}
			for (int i = 0; i < N; i++) {
				if (que.empty()) break;
				nbr.push_back(que.top().second);
				que.pop();
			}
		}

		void permutate(R range) {
			std::default_random_engine gen;
			std::uniform_real_distribution<R> dis(-range, range);
			for (int p = 0; p < np; p++) {
				if (type[p] != FLUID || surf[p] > 0.5) continue;
				pos[0][p] += dp* dis(gen);
				pos[1][p] += dp* dis(gen);
			}
		}

		void init() {
			eps = static_cast<R>(1.e-10);
			eps_mat = static_cast<R>(1.e-6);
			for (int p = 0; p < np; p++) {
				invMat.push_back(MatPP());
				invNeu.push_back(MatPP());
			}
			varrho = 1.0 / (1.0*dp);
			Vec zero = Vec::Zero();
			DR::Run<1, 0>(varrho, zero.data(), pn_px_o.data());
			DR::Run<0, 1>(varrho, zero.data(), pn_py_o.data());
			DR::Run<2, 0>(varrho, zero.data(), pn_pxx_o.data());
			DR::Run<1, 1>(varrho, zero.data(), pn_pxy_o.data());
			DR::Run<0, 2>(varrho, zero.data(), pn_pyy_o.data());
			DR::Run<1, 0>(varrho, zero.data(), pn_p_o.block<1, PN::value>(0, 0).data());
			DR::Run<0, 1>(varrho, zero.data(), pn_p_o.block<1, PN::value>(1, 0).data());
			DR::Run<2, 0>(varrho, zero.data(), pn_pp_o.block<1, PN::value>(0, 0).data());
			DR::Run<1, 1>(varrho, zero.data(), pn_pp_o.block<1, PN::value>(1, 0).data());
			DR::Run<0, 2>(varrho, zero.data(), pn_pp_o.block<1, PN::value>(2, 0).data());
			pn_lap_o = pn_pxx_o + pn_pyy_o;

			DRH::Run<1, 0>(varrho, zero.data(), pnH_px_o.data());
			DRH::Run<0, 1>(varrho, zero.data(), pnH_py_o.data());

			updateInvMat();
			makeSurf();
			
			for (int p = 0; p < np; p++) {
				div[p] = Div(vel[0].data(), vel[1].data(), p);
				vort[p] = Rot(vel[0].data(), vel[1].data(), p);
				spd[p] = sqrt(vel[0][p]*vel[0][p] + vel[1][p]*vel[1][p]);
			}
		}
		void addPart(const int& t, const Vec& p, const Vec& v) {
			Particle<R,2>::addPart(t, p, v);
			invMat.push_back(MatPP());
			invNeu.push_back(MatPP());
		}
		void erasePart(const int& offset) {
			Particle<R,2>::erasePart(offset);
			invMat.erase(invMat.begin() + offset);
			invNeu.erase(invNeu.begin() + offset);
		}

	public:
		R eps;
		R eps_mat;

		std::vector<MatPP> invMat;
		std::vector<MatPP> invNeu;

		R varrho;
		Eigen::Matrix<R, 1, PN::value, Eigen::RowMajor>						pn_px_o;
		Eigen::Matrix<R, 1, PN::value, Eigen::RowMajor>						pn_py_o;
		Eigen::Matrix<R, 1, PN::value, Eigen::RowMajor>						pn_pxx_o;
		Eigen::Matrix<R, 1, PN::value, Eigen::RowMajor>						pn_pxy_o;
		Eigen::Matrix<R, 1, PN::value, Eigen::RowMajor>						pn_pyy_o;
		Eigen::Matrix<R, 2, PN::value, Eigen::RowMajor>						pn_p_o;
		Eigen::Matrix<R, mMath::H<2, 2>::value, PN::value, Eigen::RowMajor>	pn_pp_o;
		Eigen::Matrix<R, 1, PN::value, Eigen::RowMajor>						pn_lap_o;

		Eigen::Matrix<R, 1, PNH::value, Eigen::RowMajor>					pnH_px_o;
		Eigen::Matrix<R, 1, PNH::value, Eigen::RowMajor>					pnH_py_o;
	private:
		bool python_initialized;
		bool numpy_initialized;
	};

	template <typename R, int P>
	class Particle_x<R,3,P> : public Particle<R,3> {};
