/*
* LICENCE
* copyright 2014 ~ ****
* Some rights reserved.
* Author: HUFANGYUAN
* Released under CC BY-NC
*/
//Particle.h
///defination of class Particle
///2016.4.22 fixed bug in b2normal() --- 
/// normal calculation is different when initial particles position are random
#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include "LinkCell.h"
#include "Header.h"
#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif

	template <typename R, int D>
	class Particle {};

	template <typename R>
	class Particle<R,1> {};

	template <typename R>
	class Particle<R,2> {
	public:
		typedef Eigen::Matrix<int,2,1> iVec;
		typedef Eigen::Matrix<R,2,1> Vec;
		typedef Eigen::Matrix<R,2,2> Mat;
	public:
		Particle() : ct(0), np(0) {}
		~Particle() {}

		void clean() {
			type.clear();
			pos[0].clear(); pos[1].clear();
			vel[0].clear(); vel[1].clear();
			pres.clear();
			phi.clear(); vort.clear(); div.clear();
			np = 0;
		}
		void operator >> (const std::string str) const {
			std::ofstream file(str, std::ofstream::out);
			file << std::scientific << std::setprecision(6) << ct << std::endl;
			file << std::scientific << std::setprecision(6) << dp << " " << k << std::endl;
			for (int p = 0; p < np; p++) {
				file << std::scientific << std::setprecision(6);
				file << type[p] << " " << pos[0][p] << " " << pos[1][p] << " " << vel[0][p] << " " << vel[1][p] << std::endl;
			}
			std::cout << " Writing Geo.in done. " << std::endl;
			file.close();
		}
		void operator << (const std::string str) {
			std::ifstream file(str);
			if (!file.is_open()) std::cout << " File " << str << " not found ! " << std::endl;
			file >> ct >> dp >> k;
			r0 = dp* k;
			std::string line;
			while(std::getline(file, line)) {
				std::istringstream iss(line);
				int typ;
				R pres;
				Vec p, v;
				iss >> typ;
				iss >> p[0] >> p[1];
				iss >> v[0] >> v[1];
				iss >> pres;
				addPart(int(typ), p, v, pres);
			}
			file.close();
			buildCell();
			std::cout << " Reading " << str << " done " << std::endl;
			std::cout << " Number of Partcles loaded: " << np << std::endl;
		}

		void addPart(const int& t, const Vec& p, const Vec& v, const R& ps) {
			type.push_back(t);
			pos[0].push_back(p[0]);	pos[1].push_back(p[1]);
			vel[0].push_back(v[0]); vel[1].push_back(v[1]);
			pres.push_back(ps); vort.push_back(R(0)); div.push_back(R(0)); spd.push_back(R(0));
			surf.push_back(0);
			np++;
		}
		void erasePart(const int& offset) {
			type.erase(type.begin() + offset);
			pos[0].erase(pos[0].begin() + offset); pos[1].erase(pos[1].begin() + offset);
			vel[0].erase(vel[0].begin() + offset); vel[1].erase(vel[1].begin() + offset);
			pres.erase(pres.begin() + offset); vort.erase(vort.begin() + offset); div.erase(div.begin() + offset);
			surf.erase(surf.begin() + offset);
			np--;
		}

		void buildCell() {
			BBox<R> b = BBox<R>();
			for (int p = 0; p < np; p++) {
				Vec pp;
				pp[0] = pos[0][p];
				pp[1] = pos[1][p];
				b += pp;
			}
			b.Expand(0.05);
			cell = new LinkCell<R,2>(b, r0);
			updateCell();
		}
		void updateCell() {
			cell->update(pos[0].data(), pos[1].data(), np);
		}
		void getBBox(R& left, R& right, R& bottom, R& top) const {
			cell->getBBox(left, right, bottom, top);
		}

	public:
		R ct;
		R dp, k, r0;
		int np;
		std::vector<int> type;
		std::vector<R> pos[2];
		std::vector<R> vel[2];

		std::vector<R> pres;
		std::vector<R> div;
		std::vector<R> vort;
		std::vector<R> spd;
		std::vector<R> surf;

		LinkCell<R, 2>* cell;
	};

	template <typename R>
	class Particle<R,3> {};
