/*
	Framework from here: https://github.com/SaschaWillems/Vulkan
	Some code for elephant including some shaders from here: https://prideout.net/blog/old/blog/index.html@tag=tessellation.html
	Some other shader code from: https://www.amazon.com/Graphics-Shaders-Theory-Practice-Second/dp/1568814348
	Some code from prototype
*/

#include "LatticeExample.cpp"

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
	
	// Base example
	//example = new LatticeExample(true);
	
	// Grid example
	example = new GridLatticeExample(true, 100.0f, 100.0f, 10, 10);
	
	// Cylinder Example
	//example = new CylinderLatticeExample(true, 20.0f, 50.0f, 8, 12);
	
	// Sphere Example
	//example = new SphereLatticeExample(true, 20.0f, 8, 8);

	// Non-uniform grid example
	//example = new NonUniformGridExample(true);

	// Grid where quads have angles != 90
	//example = new NonRectangularExample(true);

	// T-locus example
	//example = new TLocusExample(true);

	// Example with several lattices
	//example = new MultiLatticeExample(true);

	example->initVulkan();
	example->setupWindow(hInstance, WndProc);
	example->prepare();
	example->renderLoop();
	delete(example);
}