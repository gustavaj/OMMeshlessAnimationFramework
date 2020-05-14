/*
	Framework from here: https://github.com/SaschaWillems/Vulkan
	Some code for elephant including some shaders from here: https://prideout.net/blog/old/blog/index.html@tag=tessellation.html
	Some code from prototype
*/

#include "LatticeExample.cpp"

// TODO: Right now resources are not properly deleted when lattices are added from the menu

// TODO: Adding new local surface types are a lot of repeat work right now
// All the different pre-eval schemes are added in almost the same way but has to be done individually
// Maybe add some evaluator function that are used by all of them?

LatticeExample* example;

LRESULT CALLBACK WndProc(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
	if (example != NULL)
	{
		example->handleMessages(hWnd, uMsg, wParam, lParam);
	}
	return (DefWindowProc(hWnd, uMsg, wParam, lParam));
}

int APIENTRY WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR pCmdLine, int nCmdShow)
{
	std::cerr << "qwe" << std::endl;
	for (size_t i = 0; i < __argc; i++) { LatticeExample::args.push_back(__argv[i]); };
	
	OML::LocalSurfaceType lsType = OML::LocalSurfaceType::Quadratic_Bezier;
	OML::EvaluationMethod evalMethod = OML::EvaluationMethod::Direct;
	bool debug = true;

	// Base example
	//example = new LatticeExample(debug, lsType, evalMethod);
	
	// Grid example
	example = new GridLatticeExample(debug, 1000.0f, 1000.0f, 10, 10, lsType, evalMethod);

	// Random grid example
	//example = new RandomGridLatticeExample(debug, 100.0f, 100.0f, 10, 10, lsType, evalMethod);
	
	// Cylinder Example
	//example = new CylinderLatticeExample(debug, 20.0f, 20.0f, 1, 4, lsType, evalMethod);
	
	// Sphere Example
	//example = new SphereLatticeExample(debug, 20.0f, 8, 8, lsType, evalMethod);

	// Non-uniform grid example
	//example = new NonUniformGridExample(debug, lsType, evalMethod);

	// Grid where quads have angles != 90
	//example = new NonRectangularExample(debug, lsType, evalMethod);

	// T-locus example
	//example = new TLocusExample(debug, lsType, evalMethod);
	//example = new TLocusx4Example(debug, lsType, evalMethod);

	// Example with several lattices
	//example = new MultiLatticeExample(debug, lsType, evalMethod);

	// Example setup
	example->camera.setPerspective(60.0f, (float)example->width / (float)example->height, 1.0f, 10000.0f);
	example->camera.setRotation(glm::vec3(0.0f, 0.0f, 0.0f));
	example->camera.setTranslation(glm::vec3(0.0f, 0.0f, -500.0f));
	example->camera.movementSpeed = 100.0f;
	example->zoom = 0.0f;
	example->zoomSpeed = 100.0f;
	//example->camera.type = Camera::CameraType::firstperson;
	example->camera.type = Camera::CameraType::lookat;

	// Example settings
	example->settings.fullscreen = false;
	example->settings.overlay = true;

	OML::Timer::Start("init", "Setup example");

	example->initVulkan();
	example->setupWindow(hInstance, WndProc);
	example->prepare();
	
	OML::Timer::Stop("init", "Setup example");

	example->renderLoop();
	delete(example);
}