#include <config.h>
#include <controller.hpp>
#include <environment.hpp>
#include <tasks.hpp>
#include <utils.hpp>

#include <cmath>
#include <string>
#include <iostream>
#include <windows.h>

using namespace std;

BasicBinaryField* bbf;
Forager* fgr;
ForagingTask* ft;
NivController* nvc;

double time = 0.0;


LRESULT CALLBACK WindowProc(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
	static HDC hdcBackBuffer;
	static HBITMAP hBitmap;
	static HBITMAP hOldBitmap;
	static RECT screen;
	
	switch (uMsg)
	{
	case WM_CREATE: 
		{
			GetClientRect(hWnd, &screen);

			hdcBackBuffer = CreateCompatibleDC(NULL);
			HDC hdc = GetDC(hWnd);
			hBitmap = CreateCompatibleBitmap(hdc, screen.right, screen.bottom);
			hOldBitmap = (HBITMAP)SelectObject(hdcBackBuffer, hBitmap);
			ReleaseDC(hWnd, hdc);

			init_rng();
			bbf = new BasicBinaryField(fieldSide, neutralPct, redPct, redRwd, redPrb, bluePct, blueRwd, bluePrb, 1.0);
			fgr = new Forager(30.0, 30.0, 8.0, -1.0, 0.0, -0.2, foragerSpeed, visualConeAngle, projectionRadius, projection_dR, projection_dTheta);
			ft = new ForagingTask(bbf, fgr, (double)fieldSide, 8.0, 9.0, 100);
			nvc = new NivController();
			fgr->SetController(nvc);
		}
		break;
	case WM_PAINT:
		{
			PAINTSTRUCT ps;
			BeginPaint(hWnd, &ps);
			BitBlt(hdcBackBuffer, 0, 0, screen.right, screen.bottom, NULL, NULL, NULL, WHITENESS);

			RECT FieldView;
			FieldView.left = 10+screen.left;
			FieldView.top = 10+screen.top;
			FieldView.right = (screen.right - screen.left)/2;
			FieldView.bottom = (screen.bottom - screen.top)/2;

			RECT CameraView;
			CameraView.left = 10+screen.left;
			CameraView.top = (screen.bottom - screen.top)/2 + 10;
			CameraView.right = (screen.right - screen.left)/2;
			CameraView.bottom = screen.bottom - 10;

			Rectangle(hdcBackBuffer, FieldView.left, FieldView.top, FieldView.right, FieldView.bottom);
			Rectangle(hdcBackBuffer, CameraView.left, CameraView.top, CameraView.right, CameraView.bottom);

			char time_str[255];
			sprintf(time_str, "Time Elapsed: %.2f", time);
			TextOut(hdcBackBuffer, FieldView.right + 50, FieldView.top, time_str, strlen(time_str));

			char rwd_str[255];
			sprintf(rwd_str, "Total Reward Earned: %.2f", fgr->GetReward());
			TextOut(hdcBackBuffer, FieldView.right + 50, FieldView.top + 30, rwd_str, strlen(rwd_str));

			char trial_str[255];
			sprintf(trial_str, "Trials Complete: %d / %d", ft->CompletedTrials(), 100);
			TextOut(hdcBackBuffer, FieldView.right + 50, FieldView.top + 60, trial_str, strlen(trial_str));

			bbf->Render(hdcBackBuffer, FieldView);
			fgr->RenderCameraView(hdcBackBuffer, CameraView);
			bbf->RenderForagerEllipse(hdcBackBuffer, FieldView, fgr);
			bbf->RenderForagerShadow(hdcBackBuffer, FieldView, fgr);

			char* header = "Color Frequency Statistics:";
			char red_pct[255], blue_pct[255], gray_pct[255];
			sprintf(red_pct, "Red: %.4f", fgr->GetColorPct(RED));
			sprintf(blue_pct, "Blue: %.4f", fgr->GetColorPct(BLUE));
			sprintf(gray_pct, "Gray: %.4f", fgr->GetColorPct(GRAY));
			TextOut(hdcBackBuffer, CameraView.right + 10, CameraView.top + 10, header, strlen(header));
			TextOut(hdcBackBuffer, CameraView.right + 10, CameraView.top + 25, red_pct, strlen(red_pct));
			TextOut(hdcBackBuffer, CameraView.right + 10, CameraView.top + 40, blue_pct, strlen(blue_pct));
			TextOut(hdcBackBuffer, CameraView.right + 10, CameraView.top + 55, gray_pct, strlen(gray_pct));


			BitBlt(ps.hdc, 0, 0, screen.right, screen.bottom, hdcBackBuffer, 0, 0, SRCCOPY);
			EndPaint(hWnd, &ps);
		}
		break;
	case WM_SIZE:
		{
			screen.right = LOWORD(lParam);
			screen.bottom = HIWORD(lParam);
			SelectObject(hdcBackBuffer, hOldBitmap);
			DeleteObject(hBitmap);
			HDC hdc = GetDC(hWnd);
			hBitmap = CreateCompatibleBitmap(hdc, screen.right, screen.bottom);
			ReleaseDC(hWnd, hdc);
			SelectObject(hdcBackBuffer, hBitmap);
		}
		break;
	case WM_DESTROY:
		{
			delete ft;
			delete bbf;
			delete fgr;
			delete nvc;
			ft = NULL;
			bbf = NULL;
			fgr = NULL;
			nvc = NULL;
			free_rng();
			SelectObject(hdcBackBuffer, hOldBitmap);
			DeleteDC(hdcBackBuffer);
			DeleteObject(hBitmap);
			PostQuitMessage(0);
		}
		break;
	}
	return DefWindowProc(hWnd, uMsg, wParam, lParam);
}



int WINAPI WinMain(HINSTANCE hInstance,
				   HINSTANCE hPrevInstance,
				   LPSTR lpCmdLine,
				   int nCmdShow)
{
	WNDCLASSEX winclass;
	HWND hWnd;
	MSG msg;

	winclass.cbSize = sizeof(WNDCLASSEX);
	winclass.style = CS_HREDRAW | CS_VREDRAW;
	winclass.lpfnWndProc = WindowProc;
	winclass.cbClsExtra = 0;
	winclass.cbWndExtra = 0;
	winclass.hInstance = hInstance;
	winclass.hIcon = LoadIcon(NULL, IDI_APPLICATION);
	winclass.hCursor = LoadCursor(NULL, IDC_ARROW);
	winclass.hbrBackground = NULL;
	winclass.lpszMenuName = NULL;
	winclass.lpszClassName = g_szWindowClassName;
	winclass.hIconSm = LoadIcon(NULL, IDI_APPLICATION);
	
	if (!RegisterClassEx(&winclass))
	{
		MessageBox(NULL, "Class Registration Failed!", "Error", 0);
		return 0;
	}

	hWnd = CreateWindowEx(NULL, 
						  g_szWindowClassName, 
						  g_szApplicationName, 
						  WS_OVERLAPPEDWINDOW | WS_VISIBLE,
						  0, 
						  0, 
						  nWindow_Width, 
						  nWindow_Height, 
						  NULL, 
						  NULL, 
						  hInstance, 
						  NULL);
	if (!hWnd)
	{
		MessageBox(NULL, "Error Creating Window!", "Error", 0);
		return 0;
	}

	ShowWindow(hWnd, nCmdShow);
	UpdateWindow(hWnd);

	bool bDone = false;
	
	while (!bDone)
	{
		while (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE))
		{
			if (msg.message == WM_QUIT)
			{
				bDone = true;
			}
			else
			{
				TranslateMessage(&msg);
				DispatchMessage(&msg);
			}
		}

		//Sleep(10);
		if (ft != NULL)
			ft->Update(sim_dT);
		time += sim_dT;
		InvalidateRect(hWnd, NULL, TRUE);
		UpdateWindow(hWnd);
	}

	UnregisterClass(g_szWindowClassName, hInstance);

	
	return 0;
}