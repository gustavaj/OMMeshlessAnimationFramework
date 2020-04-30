/*
	Framework from here: https://github.com/SaschaWillems/Vulkan
	Some code for elephant including some shaders from here: https://prideout.net/blog/old/blog/index.html@tag=tessellation.html
	Some other shader code from: https://www.amazon.com/Graphics-Shaders-Theory-Practice-Second/dp/1568814348
	Some code from prototype
*/

#include "LatticeExample.cpp"

// TODO: Right now resources are not properly deleted when lattices are added from the menu

const bool debug = true;

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
	std::cerr << "asd" << std::endl;
	for (size_t i = 0; i < __argc; i++) { LatticeExample::args.push_back(__argv[i]); };
	
	OML::LocalSurfaceType lsType = OML::LocalSurfaceType::Quadratic_Bezier;
	OML::EvaluationMethod evalMethod = OML::EvaluationMethod::Pre_Sampled_Image;

	// Base example
	//example = new LatticeExample(debug, lsType, evalMethod);
	
	// Grid example
	example = new GridLatticeExample(debug, 30.0f, 30.0f, 3, 3, lsType, evalMethod);

	// Random grid example
	//example = new RandomGridLatticeExample(debug, 100.0f, 100.0f, 10, 10, lsType, evalMethod);
	
	// Cylinder Example
	//example = new CylinderLatticeExample(debug, 20.0f, 50.0f, 8, 12, lsType, evalMethod);
	
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

	OML::Timer::Start("init", "Setup example");

	example->initVulkan();
	example->setupWindow(hInstance, WndProc);
	example->prepare();
	
	OML::Timer::Stop("init", "Setup example");

	example->renderLoop();
	delete(example);
}