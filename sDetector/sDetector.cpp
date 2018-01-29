/*
* LICENCE
* copyright 2017 ~ ****
* Some rights reserved.
* Author: HUFANGYUAN
* Released under CC BY-NC
*/
// sDetector.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
///sDetector.cpp main loop
#pragma once
#include <Controller.h>
#include <VisualizationDll.h>
#include "Header.h"
#include "Bitmap.h"
#include "Particle_x.h"

struct Para {
	typedef double DataType;
	typedef DataType* DataTypePtr;
	enum { Dimension = 2, };
	enum { Order = 2, };
};

typedef void* NPtr;
typedef VIS::VisualizationDll Visualization;
Particle_x<Para::DataType, Para::Dimension, Para::Order> part;

void InterpolationWraper(const double x, const double y, double& u, double& v) {
	Eigen::Matrix<double, 2, 1> uv = part.interpolateLSA(part.vel[0].data(), part.vel[1].data(), x, y);
	u = uv[0];
	v = uv[1];
}

static VIS::Controller control;
static TwBar* GUIBar;
static TwBar* StreamBar;
static void setTwVisible(TwBar* const bar, const int visible);

static void Render() {
	switch (control.m_mode) {
	case VIS::DMODE_ZERO:
		Visualization::Run(&control, part.np, NPtr(part.type.data()), NPtr(part.pos[0].data()), NPtr(part.pos[1].data()), NPtr(part.surf.data()));
		setTwVisible(StreamBar, 0);
		break;
	case VIS::DMODE_ONE:
		Visualization::Run(&control, part.np, NPtr(part.type.data()), NPtr(part.pos[0].data()), NPtr(part.pos[1].data()), NPtr(part.vel[0].data()));
		setTwVisible(StreamBar, 0);
		break;
	case VIS::DMODE_TWO:
		Visualization::Run(&control, part.np, NPtr(part.type.data()), NPtr(part.pos[0].data()), NPtr(part.pos[1].data()), NPtr(part.vel[1].data()));
		setTwVisible(StreamBar, 0);
		break;
	case VIS::DMODE_THREE:
		Visualization::Run(&control, part.np, NPtr(part.type.data()), NPtr(part.pos[0].data()), NPtr(part.pos[1].data()), NPtr(part.spd.data()));
		setTwVisible(StreamBar, 0);
		break;
	case VIS::DMODE_FOUR:
		Visualization::Run(&control, part.np, NPtr(part.type.data()), NPtr(part.pos[0].data()), NPtr(part.pos[1].data()), NPtr(part.pres.data()));
		setTwVisible(StreamBar, 0);
		break;
	case VIS::DMODE_FIVE:
		Visualization::Run(&control, part.np, NPtr(part.type.data()), NPtr(part.pos[0].data()), NPtr(part.pos[1].data()), NPtr(part.vort.data()));
		setTwVisible(StreamBar, 0);
		break;
	case VIS::DMODE_SIX:
		Visualization::Run(&control, part.np, NPtr(part.type.data()), NPtr(part.pos[0].data()), NPtr(part.pos[1].data()), NPtr(part.div.data()));
		setTwVisible(StreamBar, 0);
		break;
	case VIS::DMODE_SEVEN:
		setTwVisible(StreamBar, 1);
		Visualization::Run_stream(&control, 0, InterpolationWraper);
		break;
	default:
		setTwVisible(StreamBar, 0);
		break;
	}
}

static void setTwVisible(TwBar* const bar, const int visible) {
	TwSetParam(bar, NULL, "visible", TW_PARAM_INT32, 1, &visible);
}

static void callBack() {
	if (control.i_print) {
		setTwVisible(GUIBar, 0);
		setTwVisible(StreamBar, 0);
		Render();
		static Bitmap bm;
		static int i = 0;
		char name[256];
		sprintf_s(name, "./out/snap%04d.png", i++);
		//bm.SaveAsBMP(name);
		bm.SaveAsPNG(name);
		setTwVisible(GUIBar, 1);
		setTwVisible(StreamBar, 1);
		control.i_print = 0;
	}
	if (control.i_sens) {
		control.i_sens = 0;
	}
	if (control.i_save) {
		static int i = 0;
		std::ostringstream oss;
		oss << "./out/" << i++ << ".surf";
		std::ofstream file(oss.str(), std::ofstream::out);
		static const int N = 8;
		for (int p = 0; p < part.np; p++) {
			if (part.type[p] == BD2) continue;
			std::vector<int> nbr;
			part.nNearestNeighbor<N>(nbr, p);
			file << std::scientific << std::setprecision(6);
			file << part.type[p] << " ";
			for (auto i : nbr) {
				file << (part.pos[0][i] - part.pos[0][p])/part.dp << " " << (part.pos[1][i] - part.pos[1][p])/part.dp << " ";
			}
			for (int i = nbr.size(); i < N; i++) file << "0 0 ";
			file << part.surf[p] << std::endl;
		}
		file.close();
		std::cout << " saved surf data. " << std::endl;
		control.i_save = 0;
	}
}
static void fps() {

}
static void onMouse(int button, int s, int x, int y) {
	if (!TwEventMouseButtonGLUT(button, s, x, y)) {
		control.clickMouse(button, s, x, y);
		if (button == GLUT_LEFT_BUTTON && s == GLUT_DOWN) {
			const int pickID = Visualization::IntersectColorPick(&control, part.np, x, y);
			if (pickID == 0x00FFFFFF) return;
			const Para::DataType* const px = (Para::DataType*)part.pos[0].data();
			const Para::DataType* const py = (Para::DataType*)part.pos[1].data();
			const Para::DataType* const div = (Para::DataType*)part.div.data();
			const Para::DataType* const pres = (Para::DataType*)part.pres.data();
			const Para::DataType* const ux = (Para::DataType*)part.vel[0].data();
			const Para::DataType* const uy = (Para::DataType*)part.vel[1].data();
			std::cout << " --------------------------------------------------------------------- " << std::endl;
			std::cout << " Particle ID : " << pickID << std::endl;
			std::cout << " Coordinate (x,y) : " << px[pickID] << ", " << py[pickID] << std::endl;
			std::cout << " Velocity (x,y) : " << ux[pickID] << ", " << uy[pickID] << std::endl;
			std::cout << " Divergence : " << div[pickID] << "    " << " Pressure : " << pres[pickID] << std::endl;
			std::cout << " --------------------------------------------------------------------- " << std::endl;
			part.surf[pickID] = part.surf[pickID] > 0.5 ? 0.0 : 1.0;
		}
	}
}
static void onMotion(int x, int y) {
	if (!TwEventMouseMotionGLUT(x, y)) {
		control.moveMouse(x, y);
		glutPostRedisplay();
	}
}
static void onMouseWheel(int button, int dir, int x, int y) {
	control.rollMouse(button, dir, x, y);
}
static void onReshape(int width, int height) {
	glViewport(0, 0, width, height);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	double left, right, bottom, top;
	part.getBBox(left, right, bottom, top);
	control.reshapeWindow(width, height, float(left), float(right), float(bottom), float(top));
	//gluPerspective(-90.0f, float(control.u_width) / float(control.u_height), 1.0f, 100.0f);
	TwWindowSize(control.u_width, control.u_height);
}
static void onKeyboard(unsigned char key, int x, int y) {
	if (!TwEventKeyboardGLUT(key, x, y)) {
		glutPostRedisplay();
		control.pressKey(key, x, y);
		///redistribute
		if (key == 'r') {
			part.Redistribute();
		}
	}
}
static void onDisplay() {
	glm::mat4 modelMatrix = glm::translate(glm::mat4(1.0f), control.m_pan)
		* glm::toMat4(control.m_rotation)
		* glm::scale(glm::mat4(1.0f), control.m_scale);

	control.m_modelMat = modelMatrix;
	control.m_viewMat = control.m_camera.GetViewMatrix();
	control.m_viewModelMat = control.m_camera.GetViewMatrix() * modelMatrix;
	control.m_projectionMat = control.m_camera.GetProjectionMatrix();
	control.m_projectionMatInv = glm::inverse(control.m_projectionMat);
	control.m_mvp = control.m_projectionMat * control.m_viewModelMat;
	control.m_mvpInv = glm::inverse(control.m_mvp);

	Render();
	TwDraw();

	glutSwapBuffers();
	glutReportErrors();

	callBack();

	glutPostRedisplay();
	if (control.i_leave) {
		glutLeaveMainLoop();
	}
}

static void onIdle() {
	return;
}

void TW_CALL ButtonSnap_callback(void*) { control.i_print = 1; }
void TW_CALL ButtonSensor_callback(void*) { control.i_sens = 1; }

void TW_CALL ButtonRenderStream_callback(void*) {
	Visualization::Run_stream(&control, 1, InterpolationWraper);
	TwDraw();
	glutSwapBuffers();
	glutReportErrors();
}

static void Run() {
	glutMainLoop();
}

static void Finalize() {
	TwTerminate();
}

static void Initialize(int argc, char** argv) {
	std::string filename;
	std::cout << " Input file name: " << std::endl;
	std::cin >> filename;
	std::cout << " Reading file... " << std::endl;
	part << filename;
	part.init();

	std::vector<Para::DataType> posx(part.pos[0]), posy(part.pos[1]);
	//part.permutate(0.3);
	for (int p = 0; p < part.np; p++) {
		posx[p] = (posx[p] - part.pos[0][p]) / part.dp;
		posy[p] = (posy[p] - part.pos[1][p]) / part.dp;
	}
	std::ofstream file("./out/" + filename, std::ofstream::out);
	static const int N = 24;
	for (int p = 0; p < part.np; p++) {
		if (part.type[p] != FLUID || part.surf[p] > 0.5) continue;
		std::vector<int> nbr;
		part.nNearestNeighbor<N>(nbr, p);
		file << std::scientific << std::setprecision(6);
		for (auto i : nbr) {
			file << (part.pos[0][i] - part.pos[0][p]) / part.dp << " " << (part.pos[1][i] - part.pos[1][p]) / part.dp << " ";
		}
		for (int i = nbr.size(); i < N; i++) {
			file << "0 0 ";
		}
		file << posx[p] << " " << posy[p] << std::endl;
	}
	std::cout << " saved redistribution vector data. " << std::endl;
	file.close();

	double left, right, bottom, top;
	part.getBBox(left, right, bottom, top);
	double whRatio = (right - left) / (top - bottom);
	control.u_height = GLuint(control.u_width / whRatio) < control.u_height_max ? GLuint(control.u_width / whRatio) : control.u_height_max;
	control.u_width = GLuint(control.u_width / whRatio) < control.u_height_max ? control.u_width : GLuint(control.u_height_max * whRatio);
	control.setProjectionOR(float(left), float(right), float(bottom), float(top));

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(control.u_width, control.u_height);
	glutCreateWindow("RTRenderer");
	glutMouseFunc(onMouse);
	glutMotionFunc(onMotion);
	glutMouseWheelFunc(onMouseWheel);
	glutReshapeFunc(onReshape);
	glutKeyboardFunc(onKeyboard);
	glutDisplayFunc(onDisplay);
	glutIdleFunc(onIdle);

	Visualization::Initialize();

	TwInit(TW_OPENGL, NULL);
	TwWindowSize(control.u_width, control.u_height);
	GUIBar = TwNewBar("GUI");
	TwDefine(" GUI size='180 300' position='0 0' ");
	TwEnumVal ev[] = { {VIS::DMODE_ZERO, "surface"}, { VIS::DMODE_ONE, "velocity-x" }, { VIS::DMODE_TWO, "velocity-y" }, { VIS::DMODE_THREE, "speed" },
					{ VIS::DMODE_FOUR, "pressure" }, { VIS::DMODE_FIVE, "vorticity" }, { VIS::DMODE_SIX, "divergence" }, {VIS::DMODE_SEVEN, "streamline"} };
	TwType quantity = TwDefineEnum("quantity", ev, 8);
	TwAddVarRW(GUIBar, "Quantity", quantity, &control.m_mode, " group='Display' ");
	TwAddVarRW(GUIBar, "Min", TW_TYPE_FLOAT, &control.f_sRangeMin, " group='Range' ");
	TwAddVarRW(GUIBar, "Max", TW_TYPE_FLOAT, &control.f_sRangeMax, " group='Range' ");
	TwDefine(" GUI/Range group='Display' ");
	TwEnumVal ev_switch[] = { { 0, "Off" }, { 1, "On" }, };
	TwAddButton(GUIBar, "Snapshot", ButtonSnap_callback, NULL, " label='Print' ");
	TwAddButton(GUIBar, "Sensor", ButtonSensor_callback, NULL, " label='Sensor' ");

	StreamBar = TwNewBar("Streamline");
	TwDefine(" Streamline size='250 200' position='180 0' ");
	TwAddVarRW(StreamBar, "Px", TW_TYPE_FLOAT, &control.v_p1.x, " group='P1' ");
	TwAddVarRW(StreamBar, "Py", TW_TYPE_FLOAT, &control.v_p1.y, " group='P1' ");
	TwAddVarRW(StreamBar, "Qx", TW_TYPE_FLOAT, &control.v_p2.x, " group='P2' ");
	TwAddVarRW(StreamBar, "Qy", TW_TYPE_FLOAT, &control.v_p2.y, " group='P2' ");
	TwAddVarRW(StreamBar, "# of streamlines", TW_TYPE_INT32, &control.i_nLines, "  ");
	TwAddVarRW(StreamBar, "Integration step size", TW_TYPE_FLOAT, &control.f_dStep, "  ");
	TwAddVarRW(StreamBar, "Streamline length", TW_TYPE_FLOAT, &control.f_sLength, "  ");
	TwAddButton(StreamBar, "Render", ButtonRenderStream_callback, NULL, " label='Render' ");
	setTwVisible(StreamBar, 0);
}

int _tmain(int argc, _TCHAR* argv[]) {
	CreateDirectoryA(std::string(".\\out").c_str(), NULL);
	Initialize(argc, (char**)argv);
	Run();
	Finalize();
	return 0;
}

