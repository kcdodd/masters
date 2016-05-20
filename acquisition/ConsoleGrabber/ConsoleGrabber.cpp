// ConsoleGrabber.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <fstream>

#define ENTER   0x0d
#define ESCAPE  0x1b 

static BOOL g_bCaptureNextFrame = FALSE;
static TCHAR g_szAppPath[_MAX_PATH + 1];

static RGBQUAD g_ColorTable[256];

static UINT g_FrameCount = 0;

void GetAppPath()
{
	_tcscpy(g_szAppPath, GetCommandLine() + 1);
	g_szAppPath[_tcsrchr(g_szAppPath, _T('\\')) - g_szAppPath] = NULL;
}

HRESULT SaveCallback(ULONG uCameraIndex, double my_double, BYTE *my_byte, BITMAPINFOHEADER *pbih)
{
    if (g_FrameCount > 0)
    {
        TCHAR       szFileName[_MAX_PATH];
        TCHAR       szSavePath[_MAX_PATH];
        SYSTEMTIME  SystemTime;

        GetLocalTime(&SystemTime);
        wsprintf(szFileName,
            _T("\\FRAME_CAM%u_%02u_%02u_%u_%02u_%02u_%02u_%03u.bmp"), uCameraIndex,
            SystemTime.wDay, SystemTime.wMonth, SystemTime.wYear, SystemTime.wHour,
            SystemTime.wMinute, SystemTime.wSecond, SystemTime.wMilliseconds);

		g_FrameCount--;

		//std::cout << SystemTime.wMinute << " " << SystemTime.wSecond << " " << SystemTime.wMilliseconds << std::endl;

		//return TRUE;

        lstrcpy(szSavePath, g_szAppPath);
        lstrcat(szSavePath, szFileName);
        //SaveBitmap(my_byte, pbih, szSavePath);
        //g_bCaptureNextFrame = FALSE;

		pbih->biCompression = BI_RGB;

		FILE *pFile = fopen(szSavePath, "wb");

		if(pFile == NULL)
        {
            cout << "Could not write bitmap" << endl;

			return TRUE;
        }
		
		BITMAPFILEHEADER bmfh;
        int nBitsOffset = sizeof(BITMAPFILEHEADER) + pbih->biSize +  256*sizeof(RGBQUAD); 
        LONG lImageSize = pbih->biSizeImage;
        LONG lFileSize = nBitsOffset + lImageSize;
        bmfh.bfType = 'B'+('M'<<8);
        bmfh.bfOffBits = nBitsOffset;
        bmfh.bfSize = lFileSize;
        bmfh.bfReserved1 = bmfh.bfReserved2 = 0;

		// *.bmp stores image upside down for some reason
		pbih->biHeight *= -1;


        //Write the bitmap file header
        UINT nWrittenFileHeaderSize = fwrite(&bmfh, 1, sizeof(BITMAPFILEHEADER), pFile);

        //And then the bitmap info header
        UINT nWrittenInfoHeaderSize = fwrite(pbih, 1, sizeof(BITMAPINFOHEADER), pFile);

		// rgb color table
		UINT nWrittenRGBQuadSize = fwrite(&g_ColorTable, 1, 256*sizeof(RGBQUAD), pFile);

		// image data
		UINT nWrittenDIBDataSize = fwrite(my_byte, 1, lImageSize, pFile);

		
		fclose(pFile);
    }

    return TRUE;
}

int _tmain(int argc, _TCHAR* argv[])
{

	int numFrames = 1;


    HRESULT hr = S_OK;
    TCHAR   szErr[128];
    CHAR    UserCmd;

	// construct color table for gray scale image
	for(size_t i = 0; i < 256; i++)
		g_ColorTable[i].rgbRed = g_ColorTable[i].rgbGreen = g_ColorTable[i].rgbBlue = (BYTE)i;

	USES_CONVERSION;

    CoInitialize(NULL);

    GetAppPath();

    CFiCamera *pMyCamera;

    pMyCamera = new CFiCamera( 0 , &hr);

    if (S_FALSE == hr)
    {
        cout << "ABORT No cameras found on the system" << endl;
		_getch();
        return -1;
    }
    else if (FAILED(hr))
    {
        wsprintf(szErr, _T("Failed to open the first camera found on the system with error %x\n"), hr);

        cout << szErr << endl;
		_getch();
        return -1;
    }

	// factory reset
	//pMyCamera->WriteQuad(pMyCamera->GetCommandRegBase() + 0x000, 1);

	// initial camera setup
	pMyCamera->SetStreamFormat(Y_MONO , res_640x480 , fps_15);
    pMyCamera->InitStream(RGB_24);
    pMyCamera->SetCallback(SaveCallback , Y_MONO);



	//pMyCamera->Set(FiExpoControl_Autoexp, 0);
	long BRIGHTNESS_OFFSET = 0x800;
	long AUTOEXP_OFFSET = 0x804;
	long SHARPNESS_OFFSET = 0x808;
	long WHITEBALANCE_OFFSET = 0x80C;
	long GAMMA_OFFSET = 0x818;
	long SHUTTER_OFFSET = 0x81C;
	long GAIN_OFFSET = 0x820;
	long IRIS_OFFSET = 0x824;
	long AM_MASK = 0x01000000;
	long cntrReg;

	long Fmt0Mode5_OFFSET = 0x214;

/*
	pMyCamera->ReadQuad(pMyCamera->GetCommandRegBase() + BRIGHTNESS_OFFSET, &cntrReg);
	cout << hex << cntrReg << endl;

	pMyCamera->ReadQuad(pMyCamera->GetCommandRegBase() + AUTOEXP_OFFSET, &cntrReg);
	cout << hex << cntrReg << endl;

	pMyCamera->ReadQuad(pMyCamera->GetCommandRegBase() + SHARPNESS_OFFSET, &cntrReg);
	cout << hex << cntrReg << endl;

	pMyCamera->ReadQuad(pMyCamera->GetCommandRegBase() + WHITEBALANCE_OFFSET, &cntrReg);
	cout << hex << cntrReg << endl;

	pMyCamera->ReadQuad(pMyCamera->GetCommandRegBase() + GAMMA_OFFSET, &cntrReg);
	cout << hex << cntrReg << endl;

	pMyCamera->ReadQuad(pMyCamera->GetCommandRegBase() + SHUTTER_OFFSET, &cntrReg);
	cout << hex << cntrReg << endl;

	pMyCamera->ReadQuad(pMyCamera->GetCommandRegBase() + GAIN_OFFSET, &cntrReg);
	cout << hex << cntrReg << endl;

	pMyCamera->ReadQuad(pMyCamera->GetCommandRegBase() + IRIS_OFFSET, &cntrReg);
	cout << hex << cntrReg << endl;
*/

	// see if autoexp is in automode
	pMyCamera->ReadQuad(pMyCamera->GetCommandRegBase() + AUTOEXP_OFFSET, &cntrReg);

	//cout << hex << cntrReg << endl;

	if (cntrReg & AM_MASK)
	{
		//cout << "turned automode off" << endl;
		// must turn off automode
		cntrReg = cntrReg & ~AM_MASK;

		pMyCamera->WriteQuad(pMyCamera->GetCommandRegBase() + AUTOEXP_OFFSET, cntrReg);
	}

	// parse camera settings file
	ifstream camsettingsfile ("camsettings.txt");
	string line;

	if (camsettingsfile.is_open())
	{
		while(!camsettingsfile.eof())
		{
			getline(camsettingsfile, line);
			int pos0 = line.find_first_not_of(" ", 0);
			int pos1 = line.find_first_of(" ", pos0);
			int pos2 = line.find_first_not_of(" ", pos1);

			string param = line.substr(pos0, pos1 - pos0);
			string svalue = line.substr(pos2, line.length() - pos2);
		
			float fvalue = atof(svalue.c_str());
			float min, max;

			if (param.compare("SHUTTER") == 0)
			{
				

				// see if it is in automode
				pMyCamera->ReadQuad(pMyCamera->GetCommandRegBase() + SHUTTER_OFFSET, &cntrReg);

				if (cntrReg & AM_MASK)
				{
					//cout << "turned automode off" << endl;
					// must turn off automode
					cntrReg = cntrReg & ~AM_MASK;

					pMyCamera->WriteQuad(pMyCamera->GetCommandRegBase() + SHUTTER_OFFSET, cntrReg);
				}

				// set shutter (rel:0 - 511. abs: )
				pMyCamera->GetRange(FiExpoControl_Shutter, &min, &max);
				std::cout << "shutter range: " << min << " , " << max << "; ";

				hr = pMyCamera->Set(FiExpoControl_Shutter, fvalue); 

				if (hr == S_FALSE || FAILED(hr))
				{
					std::cout << "Could not set shutter" << std::endl;
				}else{
					pMyCamera->Get(FiExpoControl_Shutter, &fvalue);
					std::cout << "Set shutter to: " << fvalue << std::endl;
				}
			}else if (param.compare("GAIN") == 0) {

				// see if it is in automode
				pMyCamera->ReadQuad(pMyCamera->GetCommandRegBase() + GAIN_OFFSET, &cntrReg);

				if (cntrReg & AM_MASK)
				{
					//cout << "turned automode off" << endl;
					// must turn off automode
					cntrReg = cntrReg & ~AM_MASK;

					pMyCamera->WriteQuad(pMyCamera->GetCommandRegBase() + GAIN_OFFSET, cntrReg);
				}

				// set gain (0 - 255)
				pMyCamera->GetRange(FiExpoControl_Gain, &min, &max);

				std::cout << "gain range: " << min << " , " << max << "; ";

				hr = pMyCamera->Set(FiExpoControl_Gain, fvalue); 

				if (hr == S_FALSE || FAILED(hr))
				{
					std::cout << "Could not set gain" << std::endl;
				}else{
					pMyCamera->Get(FiExpoControl_Gain, &fvalue);
					std::cout << "Set gain to: " << fvalue << std::endl;
				}
			}else if (param.compare("IRIS") == 0) {

				// see if it is in automode
				pMyCamera->ReadQuad(pMyCamera->GetCommandRegBase() + IRIS_OFFSET, &cntrReg);

				if (cntrReg & AM_MASK)
				{
					//cout << "turned automode off" << endl;
					// must turn off automode
					cntrReg = cntrReg & ~AM_MASK;

					pMyCamera->WriteQuad(pMyCamera->GetCommandRegBase() + IRIS_OFFSET, cntrReg);
				}

				// set iris 
				pMyCamera->GetRange(FiExpoControl_Iris, &min, &max);

				std::cout << "iris range: " << min << " , " << max << "; ";

				hr = pMyCamera->Set(FiExpoControl_Iris, fvalue); 

				if (hr == S_FALSE || FAILED(hr))
				{
					std::cout << "Could not set iris" << std::endl;
				}else{
					pMyCamera->Get(FiExpoControl_Iris, &fvalue);
					std::cout << "Set iris to: " << fvalue << std::endl;
				}

			}else if (param.compare("BRIGHTNESS") == 0) {

				// see if it is in automode
				pMyCamera->ReadQuad(pMyCamera->GetCommandRegBase() + BRIGHTNESS_OFFSET, &cntrReg);

				if (cntrReg & AM_MASK)
				{
					//cout << "turned automode off" << endl;
					// must turn off automode
					cntrReg = cntrReg & ~AM_MASK;

					pMyCamera->WriteQuad(pMyCamera->GetCommandRegBase() + BRIGHTNESS_OFFSET, cntrReg);
				}

				// set brightness (128 - 383)
				pMyCamera->GetRange(FiBasicControl_Brightness, &min, &max);

				std::cout << "brightness range: " << min << " , " << max << "; ";

				hr = pMyCamera->Set(FiBasicControl_Brightness, fvalue); 

				if (hr == S_FALSE || FAILED(hr))
				{
					std::cout << "Could not set brightness" << std::endl;
				}else{
					pMyCamera->Get(FiBasicControl_Brightness, &fvalue);
					std::cout << "Set brightness to: " << fvalue << std::endl;
				}

			}else if (param.compare("SHARPNESS") == 0) {

				// see if it is in automode
				pMyCamera->ReadQuad(pMyCamera->GetCommandRegBase() + SHARPNESS_OFFSET, &cntrReg);

				if (cntrReg & AM_MASK)
				{
					//cout << "turned automode off" << endl;
					// must turn off automode
					cntrReg = cntrReg & ~AM_MASK;

					pMyCamera->WriteQuad(pMyCamera->GetCommandRegBase() + SHARPNESS_OFFSET, cntrReg);
				}

				// set sharpness (0, 255)
				pMyCamera->GetRange(FiBasicControl_Sharpness, &min, &max);

				std::cout << "sharpness range: " << min << " , " << max << "; ";


				hr = pMyCamera->Set(FiBasicControl_Sharpness, fvalue); 

				if (hr == S_FALSE || FAILED(hr))
				{
					std::cout << "Could not set sharpness" << std::endl;
				}else{
					pMyCamera->Get(FiBasicControl_Sharpness, &fvalue);
					std::cout << "Set sharpness to: " << fvalue << std::endl;
				}
			}else if (param.compare("GAMMA") == 0) {

				// see if it is in automode
				pMyCamera->ReadQuad(pMyCamera->GetCommandRegBase() + GAMMA_OFFSET, &cntrReg);

				if (cntrReg & AM_MASK)
				{
					//cout << "turned automode off" << endl;
					// must turn off automode
					cntrReg = cntrReg & ~AM_MASK;

					pMyCamera->WriteQuad(pMyCamera->GetCommandRegBase() + GAMMA_OFFSET, cntrReg);
				}

				// set gamma (0 - 1)
				pMyCamera->GetRange(FiBasicControl_Gamma, &min, &max);

				std::cout << "gamma range: " << min << " , " << max << "; ";

				hr = pMyCamera->Set(FiBasicControl_Gamma, fvalue); 

				if (hr == S_FALSE || FAILED(hr))
				{
					std::cout << "Could not set gamma" << std::endl;
				}else{
					pMyCamera->Get(FiBasicControl_Gamma, &fvalue);
					std::cout << "Set gamma to: " << fvalue << std::endl;
				}
			}else if (param.compare("FPS") == 0) {
				bool ok = false;
				std::cout << "Accepted FPS: ";

				pMyCamera->ReadQuad(pMyCamera->GetCommandRegBase() + Fmt0Mode5_OFFSET, &cntrReg);

				if (cntrReg & 0x40000000){
					std::cout << "3.75, ";

					if (fvalue == 3.75)
						ok = true;
				}

				if (cntrReg & 0x20000000){
					std::cout << "7.5, ";

					if (fvalue == 7.5)
						ok = true;
				}

				if (cntrReg & 0x10000000){
					std::cout << "15, ";

					if (fvalue == 15)
						ok = true;
				}

				if (cntrReg & 0x08000000){
					std::cout << "30, ";

					if (fvalue == 30)
						ok = true;
				}

				if (cntrReg & 0x04000000){
					std::cout << "60, ";

					if (fvalue == 60)
						ok = true;
				}

				std::cout << "; ";

				// set frames per second
				if (!ok){
					std::cout << "Could not set fps" << std::endl;
				}else if (fvalue == 3.75){
					pMyCamera->SetStreamFormat(Y_MONO , res_640x480 , fps_3_75);
					std::cout << "Set fps to: 3.75" << std::endl;
				}else if (fvalue == 7.5){
					pMyCamera->SetStreamFormat(Y_MONO , res_640x480 , fps_7_5);
					std::cout << "Set fps to: 7.5" << std::endl;
				}else if (fvalue == 15){
					pMyCamera->SetStreamFormat(Y_MONO , res_640x480 , fps_15);
					std::cout << "Set fps to: 15" << std::endl;
				}else if (fvalue == 30){
					pMyCamera->SetStreamFormat(Y_MONO , res_640x480 , fps_30);
					std::cout << "Set fps to: 30" << std::endl;
				}else if (fvalue == 60){
					pMyCamera->SetStreamFormat(Y_MONO , res_640x480 , fps_60);
					std::cout << "Set fps to: 60" << std::endl;
				}
			}else if (param.compare("NUMFRAMES") == 0) {
				// set number of frames
				numFrames = fvalue;
				std::cout << "Set num-frames to: " << numFrames << std::endl;
			}

		}

		camsettingsfile.close();

	}

	


    pMyCamera->Run();

	_tprintf(_T("%s %s %u Camera started\n"), 
		A2T(pMyCamera->m_VendorInfo.szCameraVendor) , A2T(pMyCamera->m_VendorInfo.szCameraModelName),
		pMyCamera->m_VendorInfo.uCameraSerial);

	FIREi_FPS fps;

	pMyCamera->GetCurrentFrameRate(&fps);

	if (fps == fps_15)
		std::cout << "fps is 15" << std::endl;

	// start user trigger input
    for(;;)
    {
        cout << "ENTER  Take Snapshot" << endl;
        cout << "ESC    Exit" << endl;
        cout <<  endl;

        UserCmd = _getch();

        switch (UserCmd)
        {
            case ENTER:
				_tprintf(_T("Saving camera snapshot at %s\n\n"), g_szAppPath);
                //g_bCaptureNextFrame = TRUE;
				g_FrameCount = numFrames;
                break;
            //------------------
            case ESCAPE:
                pMyCamera->Stop();
                pMyCamera->ShutdownStream();
                delete pMyCamera;
                CoUninitialize( );
                return 0;
            //------------------
            default:
                cout << "Unsupported Command" << endl;
                break;
        }
    }

    pMyCamera->Stop();
    pMyCamera->ShutdownStream();

    delete pMyCamera;

    CoUninitialize( );
}

