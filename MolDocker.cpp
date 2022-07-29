
#include "framework.h"

float dt = 0.01f;

//GRAFIKA VIDEO 5.5
class Camera2D {
private:
	vec2 wCenter;
	vec2 wSize;
public:
	Camera2D() : wCenter(0, 0), wSize(600, 600) {}
	mat4 View() { return TranslateMatrix(-wCenter); }
	mat4 Project() { return ScaleMatrix(vec2(2 / wSize.x, 2 / wSize.y)); }

	mat4 ViewInv() { return TranslateMatrix(wCenter); }
	mat4 ProjectInv() { return ScaleMatrix(vec2(wSize.x / 2, wSize.y / 2)); }

	void zoom(float s) { wSize = wSize * s; }
	void pan(vec2 t) { wCenter = wCenter + t; }
};

Camera2D cam;
GPUProgram gpuProgram;

class Circle {
private:
	unsigned int vao;
	int nv;
	float phi = 0.0f;
	float sx, sy;
	float r = 0.0f;
	vec2 cgrav;
	vec2 wTranslate;
	vec3 color;
public:
	Circle(float r, vec3 c) : sx(0), sy(0), r(r), color(c), nv(100) {}
	void setTrans(vec2 wT) {
		wTranslate = wTranslate + wT;
	}
	void create() {
		const int nv = 150;
		this->nv = nv;
		glGenVertexArrays(1, &vao); //Vertex array object creation
		glBindVertexArray(vao);		//Activating it

		unsigned int vbo;
		glGenBuffers(1, &vbo);
		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		vec2 vertices[nv];
		for (int i = 0; i < nv; i++) {
			float fi = i * 2 * M_PI / nv;
			vertices[i] = vec2(cosf(fi), sinf(fi));
		}
		glBufferData(GL_ARRAY_BUFFER, sizeof(vec2) * nv, vertices, GL_DYNAMIC_DRAW);

		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, NULL);
	}
	mat4 M() {
		mat4 mScale(r + sx, 0, 0, 0,
			0, r + sy, 0, 0,
			0, 0, 0, 0,
			0, 0, 0, 1);

		mat4 trans(1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			cgrav.x, cgrav.y, 0, 1);

		mat4 rot(cosf(phi), sinf(phi), 0, 0,
			-sinf(phi), cosf(phi), 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1);

		mat4 beforeTrans(1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			wTranslate.x - cgrav.x, wTranslate.y - cgrav.y, 0, 1);

		return mScale * beforeTrans * rot * trans;
	}
	void draw() {
		int loc = glGetUniformLocation(gpuProgram.getId(), "color");
		glUniform3f(loc, color.x, color.y, color.z);
		mat4 MVPtransform = M() * cam.View() * cam.Project();
		gpuProgram.setUniform(MVPtransform, "MVP");

		glBindVertexArray(vao);
		glDrawArrays(GL_TRIANGLE_FAN, 0, nv);
	}

	vec2 getWTrans() { 
		return wTranslate; 
	}
	void addRotation(float angle) { phi += angle; }
	void setColor(vec3 c) { color = c; }
	void setCGrav(vec2 g) { cgrav = g; }
};

class StrippedLine {
private:
	float phi = 0.0f;
	unsigned int vao, vbo;
	vec2 wTranslate;
	vec2 cgrav;
	std::vector<float> vertexData;
public:
	void create(vec2 from, vec2 to) {
		glGenVertexArrays(1, &vao);
		glBindVertexArray(vao);

		glGenBuffers(1, &vbo);
		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		
		int numberOfSections = 100;
		for (int i = 0; i < numberOfSections; i++) {
			vertexData.push_back(from.x + i * (to.x - from.x) / numberOfSections);
			vertexData.push_back(from.y + i * (to.y - from.y) / numberOfSections);
		}
		glBufferData(GL_ARRAY_BUFFER, vertexData.size() * sizeof(float), &vertexData[0], GL_DYNAMIC_DRAW);

		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), reinterpret_cast<void*>(0));
	}

	mat4 M() {
		mat4 trans(1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			cgrav.x, cgrav.y, 0, 1);

		mat4 beforeTrans(1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			wTranslate.x - cgrav.x, wTranslate.y - cgrav.y, 0, 1);

		mat4 rot(cosf(phi), sinf(phi), 0, 0,
			-sinf(phi), cosf(phi), 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1);

		return beforeTrans * rot * trans;
	}

	void addTranslation(vec2 wt) {
		wTranslate = wTranslate + wt;
	}
	void addRotation(float angle) { phi += angle; }
	void setCGrav(vec2 g) { cgrav = g; }
	void draw() {
		if (vertexData.size() > 0) {
			int loc = glGetUniformLocation(gpuProgram.getId(), "color");
			glUniform3f(loc, 1.0f, 1.0f, 1.0f);
			mat4 MVPtransform = M() * cam.View() * cam.Project();
			gpuProgram.setUniform(MVPtransform, "MVP");
			glBindVertexArray(vao);
			glDrawArrays(GL_LINE_STRIP, 0, vertexData.size() / 2);
		}
	}
};

class Atom {
private:
	Circle* circle;
	vec2 force;
	StrippedLine* line;
	float charge = 0.0f;
	float weight = 0.0f;
	float rad	 = 0.0f;

	void calculateColor() {
		double visualCharge = charge * 0.1;
		if (charge < 0) {
			visualCharge *= -1;
			if (visualCharge > 1.0f)
				visualCharge = 0.9f;
			circle->setColor(vec3(0.0f, 0.0f,0.1 + visualCharge));
		}
		if (charge > 0) {
			if (visualCharge > 1.0f)
				visualCharge = 0.9f;
			circle->setColor(vec3(0.1f + visualCharge, 0.0f, 0.0f));
		}
	}

public:
	Atom(double c, double w, double r, vec2 cen): charge(c), weight(w), rad(r) {
		circle = new Circle(rad, vec3(0, 0, 0));
		line = new StrippedLine();
		calculateColor();
		circle->create();
		circle->setTrans(cen);
	}
	~Atom() {
		delete line;
		delete circle;
	}
	void drawAtom() {
		circle->draw();
	}
	void drawConnection() {
		line->draw();
	}

	float getCharge() { return charge; }
	float getWeight() { return weight; }
	vec2 getCenter() {
		return circle->getWTrans();
	}
	vec2 getForce() { return force; }
	Circle* getCircle() { return circle; }
	StrippedLine* getLine() { return line; }

	void setCenter(vec2 wt) { circle->setTrans(wt); }
	void setForce(vec2 f) { force = f; }


	void setConnection(Atom* neigh) {
		line->create(getCenter(), neigh->getCenter());
	}
	bool doesInterlope(vec2 other, float rad) {
		return this->rad + rad > (sqrt(pow(getCenter().x - other.x, 2) + pow(getCenter().y - other.y, 2)));
	}
};

class Molecule {
private:
	Atom* atoms[8];
	int numberOfAtoms;
	float torque = 0;
	float inertia = 0;
	float angular = 0;
	float weightSum = 0;
	vec2 centerOfGravity;
	vec2 forceSum;
	vec2 speed;
	vec2 startingPoint;

	float randomWeight() {
		return (std::rand() % 5 + 1) * 1.66e-8; //nagyságrendek tartása
	}
	float randomCharge(bool pos) {
		return (std::rand() % 5 + 1)* 1.6 * (pos ? 1 : -1);
	}
	float randomRad(float weight) {
		return weight * 1.66e8 * 2;
	}
	vec2 randomLoc(float rad, int atomNum) {
		vec2 re = vec2(std::rand() % 300 + startingPoint.x, std::rand() % 300 + startingPoint.y);
		while (true) {
			bool goodPosition = true;
			for (int i = 0; i < atomNum; i++) {
				if (atoms[i]->doesInterlope(re, rad)) {
					goodPosition = false;
					break;
				}
			}
			if (goodPosition)
				return re;
			re = vec2(std::rand() % 300 + startingPoint.x, std::rand() % 300 + startingPoint.y);
		}

	}

	void genererateAtoms() {
		numberOfAtoms = std::rand() % 7 + 2;
		float totalCharge = 0;
		for (int i = 0; i < numberOfAtoms; i++) {
			float charge = randomCharge((std::rand() % 2) ? true : false);
			float weight = randomWeight();
			float rad = randomRad(weight);
			if (lastAtom(i)) {
				charge = totalCharge * -1;
			}
			else
				totalCharge += charge;
			atoms[i] = new Atom(charge, weight, rad, randomLoc(rad, i));
			if (i > 0) {
				atoms[i]->setConnection(atoms[std::rand() % i]);
			}
			weightSum += atoms[i]->getWeight();
		}
		
	}
	void calculateCenterOfGravity() {
		centerOfGravity = vec2(0, 0);
		vec2 refPoint = vec2(atoms[0]->getCenter());
		double totalWeight = 0;
		for (int i = 0; i < numberOfAtoms; i++) {
			totalWeight += atoms[i]->getWeight();
			centerOfGravity = centerOfGravity + atoms[i]->getCenter() * atoms[i]->getWeight();
		}
		centerOfGravity = centerOfGravity / totalWeight;
		for (int i = 0; i < numberOfAtoms; i++) {
			atoms[i]->getCircle()->setCGrav(centerOfGravity);
			atoms[i]->getLine()->setCGrav(centerOfGravity);
		}
	}

	/* Conditions */
	bool lastAtom(int i) { return i == numberOfAtoms -1; }
public:
	Molecule(vec2 startingPoint) : startingPoint(startingPoint) {
		genererateAtoms();

		calculateCenterOfGravity();
	}
	~Molecule() {
		for (int i = 0; i < numberOfAtoms; i++)
			delete atoms[i];
	}
	void draw() {
		for (int i = 0; i < numberOfAtoms; i++)
			atoms[i]->drawConnection();
		for (int i = 0; i < numberOfAtoms; i++)
			atoms[i]->drawAtom();
	}
	void moveMolecule(Molecule* other) {
		calcForce(other);
		calcSpeed();
		moveAtoms();
		calculateCenterOfGravity();
		calcTorque();
		calcInertia();
		rotateAtoms();
		glutPostRedisplay();
	}
	void calcForce(Molecule* other) {
		forceSum = vec2(0, 0);
		for (int i = 0; i < numberOfAtoms; i++) {
			atoms[i]->setForce(0);
			for (int j = 0; j < other->getAtomNumber(); j++) {
				vec2 fromTo = vec2(atoms[i]->getCenter().x - other->getAtoms()[j]->getCenter().x, atoms[i]->getCenter().y - other->getAtoms()[j]->getCenter().y);
				if (fromTo.x != 0.0 && fromTo.y != 0.0) {
					atoms[i]->setForce(atoms[i]->getForce() + (atoms[i]->getCharge() * other->getAtoms()[j]->getCharge() / (2 * M_PI * 8.854e7 * length(fromTo)) * normalize(fromTo)) * 10e9);
				}
			}
			forceSum = forceSum + atoms[i]->getForce();
		}
		forceSum = forceSum - speed * 3;
	}
	void calcTorque() {
		torque = 0;
		for (int i = 0; i < numberOfAtoms; i++) {
			vec2 fromTo = atoms[i]->getCenter() - centerOfGravity;
			vec3 z = cross(fromTo, atoms[i]->getForce());
			torque = torque + z.z;
		}
	}
	void calcSpeed() {
		speed = speed + ((forceSum / weightSum) * dt / 10e5);
	}
	void calcInertia() {
		inertia = 0;
		for (int i = 0; i < numberOfAtoms; i++) {
			vec2 fromTo = atoms[i]->getCenter() - centerOfGravity;
			inertia += atoms[i]->getWeight() * pow(length(fromTo), 2);
		}
	}
	void moveAtoms() {
		for (int i = 0; i < numberOfAtoms; i++) {
			atoms[i]->setCenter(speed * dt);
			atoms[i]->getLine()->addTranslation(speed * dt);
		}
	}
	void rotateAtoms() {
		angular += torque / inertia * dt;
		float angle = angular * dt / 9.2e8;
		for (int i = 0; i < numberOfAtoms; i++) {
			atoms[i]->getLine()->addRotation(angle);
			atoms[i]->getCircle()->addRotation(angle);
		}
		//l.phi += 0.3;
	}

	int getAtomNumber() { return numberOfAtoms; }
	Atom** getAtoms() { return atoms; }
	vec2 getSpeed() { return speed; }
};

Molecule* m1 = nullptr;
Molecule* m2 = nullptr;


// vertex shader in GLSL: It is a Raw string (C++11) since it contains new line characters
const char * const vertexSource = R"(
	#version 330				// Shader 3.3
	precision highp float;		// normal floats, makes no difference on desktop computers

	uniform mat4 MVP;			// uniform variable, the Model-View-Projection transformation matrix
	layout(location = 0) in vec2 vp;	// Varying input: vp = vertex position is expected in attrib array 0

	void main() {
		vec4 re = vec4(vp.x, vp.y, 0, 1) * MVP;
		float t = 1 / (sqrt(re.x * re.x + re.y * re.y + 1) + 1);
		re.x *= t;
		re.y *= t;
		gl_Position = vec4(re.x, re.y, 0, 1);		// transform vp from modeling space to normalized device space
	}
)";

// fragment shader in GLSL
const char * const fragmentSource = R"(
	#version 330			// Shader 3.3
	precision highp float;	// normal floats, makes no difference on desktop computers
	
	uniform vec3 color;		// uniform variable, the color of the primitive
	out vec4 outColor;		// computed color of the current pixel

	void main() {
		outColor = vec4(color, 1);	// computed color is the color of the primitive
	}
)";


// Initialization, create an OpenGL context
void onInitialization() {
	glViewport(0, 0, windowWidth, windowHeight);
	glLineWidth(2.0f);

	m1 = new Molecule(vec2(std::rand() % 500 - 250, 0.0));
	m2 = new Molecule(vec2(std::rand() % 500 - 250, 0.0));


	gpuProgram.create(vertexSource, fragmentSource, "outColor");
}

// Window has become invalid: Redraw
void onDisplay() {
	glClearColor(0.2, 0.2, 0.2, 1);     // background color
	glClear(GL_COLOR_BUFFER_BIT); // clear frame buffer

	m1->draw();
	m2->draw();


	glutSwapBuffers(); // exchange buffers for double buffering
}

void spawnTwoMolecule();
// Key of ASCII code pressed
void onKeyboard(unsigned char key, int pX, int pY) {
	switch (key) {
	case ' ': spawnTwoMolecule(); break;
	case 'z': cam.zoom(0.9f); break;
	case 'o': cam.zoom(1.1f); break;
	case 'd': cam.pan(vec2( 1,  0)); break;
	case 'e': cam.pan(vec2( 0,  1)); break;
	case 's': cam.pan(vec2(-1,  0)); break;
	case 'x': cam.pan(vec2( 0, -1)); break;
	}
	glutPostRedisplay();
}

// Key of ASCII code released
void onKeyboardUp(unsigned char key, int pX, int pY) {
}

// Move mouse with key pressed
void onMouseMotion(int pX, int pY) {	// pX, pY are the pixel coordinates of the cursor in the coordinate system of the operation system
	// Convert to normalized device space
	float cX = 2.0f * pX / windowWidth - 1;	// flip y axis
	float cY = 1.0f - 2.0f * pY / windowHeight;
	printf("Mouse moved to (%3.2f, %3.2f)\n", cX, cY);
}

// Mouse click event
void onMouse(int button, int state, int pX, int pY) { // pX, pY are the pixel coordinates of the cursor in the coordinate system of the operation system
	// Convert to normalized device space
	float cX = 2.0f * pX / windowWidth - 1;	// flip y axis
	float cY = 1.0f - 2.0f * pY / windowHeight;

	char * buttonStat;
	switch (state) {
	case GLUT_DOWN: buttonStat = "pressed"; break;
	case GLUT_UP:   buttonStat = "released"; break;
	}

	switch (button) {
	case GLUT_LEFT_BUTTON:   printf("Left button %s at (%3.2f, %3.2f)\n", buttonStat, cX, cY);   break;
	case GLUT_MIDDLE_BUTTON: printf("Middle button %s at (%3.2f, %3.2f)\n", buttonStat, cX, cY); break;
	case GLUT_RIGHT_BUTTON:  printf("Right button %s at (%3.2f, %3.2f)\n", buttonStat, cX, cY);  break;
	}
}
// Idle event indicating that some time elapsed: do animation here
void onIdle() {
	long time = glutGet(GLUT_ELAPSED_TIME); // elapsed time since the start of the program
	if (m1 != nullptr && m2 != nullptr) {
		m1->moveMolecule(m2);
		m2->moveMolecule(m1);
	}
}

void spawnTwoMolecule() {
	delete m1;
	delete m2;
	m1 = new Molecule(vec2(std::rand() % 100-50, 0.0));
	m2 = new Molecule(vec2(std::rand() % 100-50, 0.0));
}
