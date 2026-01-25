//mingw32-make -j (bedąc w build)

// Funkcja dyskretyzowana -- Przycisk h
// Wczytywanie pliku mp3 z folderu -- Przycisk m
// Po naciśnięciu m należy w terminalu wybrac przyskis do
// przeanalizowania porzez program w terminalu.

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <fstream>
#include <cstdint>

#include "CLI/CLI.hpp"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include "plotter.h"
#include "overlay.h"

#include <complex>
#include <valarray>
#define DR_MP3_IMPLEMENTATION
#include "dr_mp3.h"
#include <filesystem>

namespace fs = std::filesystem;

const double PI = 3.14159265358979323846;
typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

/////////////////////////////////////////////////////////////////////////

PlotData plot_data;

/////////////////////////////////////////////////////////////////////////

namespace {

const char* vertexShaderSource = R"(
#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec2 aTexCoord;
out vec2 TexCoord;
uniform mat4 projection;
void main()
{
    gl_Position = projection * vec4(aPos, 1.0);
    TexCoord = aTexCoord;
}
)";
const char* fragmentShaderSource = R"(
#version 330 core
out vec4 FragColor;
in vec2 TexCoord;
uniform sampler2D ourTexture;
void main()
{
    FragColor = texture(ourTexture, TexCoord);
}
)";

const char* vertexOvlSource = R"(
#version 330 core
layout (location = 0) in vec3 aPos;
uniform mat4 projection;
void main()
{
    gl_Position = projection * vec4(aPos, 1.0);
}
)";
const char* fragmentOvlSource = R"(
#version 330 core
out vec4 FragColor;
uniform vec3 objectColor;
void main()
{
    FragColor = vec4(objectColor, 1.0);
}
)";

int texture_width_ = 0;
int texture_height_ = 0;
unsigned int texture_id_;

std::unique_ptr<RenderObject> points_;
std::unique_ptr<RenderObject> lines_x_;
std::unique_ptr<RenderObject> lines_y_;

int key_pressed_ = 0;

const std::string plot_filename_ = "plot.png";

} // end of anonymous namespace

/////////////////////////////////////////////////////////////////////////

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void key_callback(GLFWwindow* window, int key, int scancode, int action,
                  int mods);
void mouse_button_callback(GLFWwindow* window, int button, int action,
                           int mods);
void cursor_position_callback(GLFWwindow* window, double xpos, double ypos);
void reloadTexture(GLFWwindow* window, const std::string& filename);

/////////////////////////////////////////////////////////////////////////

// plot to norm plot
Vertex plot_to_normal_coords(const Vertex & v)
{
    return Vertex{ (v.x - plot_data.range_x_min)
                       / (plot_data.range_x_max - plot_data.range_x_min),
                   (v.y - plot_data.range_y_min)
                       / (plot_data.range_y_max - plot_data.range_y_min) };
}

// norm plot to pixel
Vertex normal_to_pixel_coords(const Vertex& n)
{
    return Vertex{
        float(plot_data.pix_x - plot_data.pad_x * 2) * n.x + plot_data.pad_x,
        float(plot_data.pix_y - plot_data.pad_y * 2) * n.y + plot_data.pad_y
    };
}

// pix to ortho
Vertex pixel_to_ortho_coords(const Vertex& v)
{
    return Vertex{ float(v.x / plot_data.pix_x) * 2.f - 1.f,
                   float(v.y / plot_data.pix_y) * 2.f - 1.f };
}

// plot to ortho coord
Vertex plot_to_ortho_coords(const Vertex & v)
{
    Vertex w = plot_to_normal_coords(v);
    w = normal_to_pixel_coords(w);
    return pixel_to_ortho_coords(w);
}

Vertex pixel_to_plot(const Vertex & v)
{
    return Vertex{ plot_data.range_x_min
                       + float(v.x - plot_data.pad_x)
                           / (plot_data.pix_x - plot_data.pad_x * 2)
                           * (plot_data.range_x_max - plot_data.range_x_min),
                   plot_data.range_y_min
                       + (1.f
                          - (float(v.y - plot_data.pad_y)
                             / (plot_data.pix_y - plot_data.pad_y * 2)))
                           * (plot_data.range_y_max - plot_data.range_y_min) };
}

/////////////////////////////////////////////////////////////////////////
// Example of creating plot from given function

double f(double x) {
    return x*x*x + x*x - 3*x - 3;
}

// Pochodna funkcji f'(x)
double df(double x) {
    return 3*x*x + 2*x - 3;
}

double secant(double x0, double x1, double tol, int max_iter) {
    double x2;
    int iter = 0;
    do {
        double f0 = f(x0);
        double f1 = f(x1);
        if (f1 - f0 == 0) {
            return x1;
        }
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
        if (fabs(x2 - x1) < tol) break;
        x0 = x1;
        x1 = x2;
        iter++;
    } while (iter < max_iter);
    return x2;
}

std::function<double(double)> GenerateLagrangePolynomialFromPoints(
    const std::unique_ptr<RenderObject>& points_)
{
    std::vector<float> X, Y;
    X.reserve(points_->vertices.size());
    Y.reserve(points_->vertices.size());

    for (const auto& v : points_->vertices)
    {
        X.push_back(v.x);
        Y.push_back(v.y);
    }

    int n = X.size();
    if (n < 2)
    {
        std::cout << "Need at least 2 points for interpolation.\n";
        return [](double){ return 0.0; };
    }

    auto addPoly = [&](const std::vector<double>& A,
                       const std::vector<double>& B)
    {
        int m = std::max(A.size(), B.size());
        std::vector<double> C(m, 0.0);

        for (int i = 0; i < A.size(); i++) C[i] += A[i];
        for (int i = 0; i < B.size(); i++) C[i] += B[i];
        return C;
    };

    auto mulByLinear = [&](const std::vector<double>& A, double a)
    {
        int m = A.size();
        std::vector<double> C(m + 1, 0.0);

        for (int i = 0; i < m; i++) {
            C[i]     -= A[i] * a;
            C[i + 1] += A[i];
        }
        return C;
    };

    auto mulByConst = [&](const std::vector<double>& A, double c)
    {
        std::vector<double> C = A;
        for (auto &v : C) v *= c;
        return C;
    };

    //obliczenie wielomianu Lagrange’a

    std::vector<double> P(n, 0.0);

    for (int i = 0; i < n; i++)
    {
        std::vector<double> Li = {1.0};
        double denom = 1.0;

        for (int j = 0; j < n; j++) {
            if (i == j) continue;
            Li = mulByLinear(Li, X[j]);
            denom *= (X[i] - X[j]);
        }

        Li = mulByConst(Li, Y[i] / denom);
        P = addPoly(P, Li);
    }

    std::cout << "P(x) = ";
    for (int i = 0; i < P.size(); i++)
        std::cout << std::showpos << P[i] << "*x^" << i << " ";
    std::cout << std::endl;


    return [P](double x) {
        double sum = 0.0;
        double xp = 1.0;
        for (double c : P) {
            sum += c * xp;
            xp *= x;
        }
        return sum;
    };
}


void on_key_z_pressed(GLFWwindow* window) 
{

    plot_data.rgb[0] = 0.0;
    plot_data.rgb[1] = 0.0;
    plot_data.rgb[2] = 1.0;

    auto poly = GenerateLagrangePolynomialFromPoints(points_);

    GeneratePlotFromFunc("plot.png", poly, 64, 1, 3);


    reloadTexture(window, "plot.png");

    updateRenderObject(points_.get());

}

void on_key_a_pressed(GLFWwindow* window)
{
    double x0_secant = 1.0; 
    double x1_secant = 2.0;  
    double precision = 1e-6;
    int max_iter = 100;

    float px = (float)secant(x0_secant, x1_secant, precision, max_iter);

    plot_data.rgb[0] = 0.0;
    plot_data.rgb[1] = 0.0;
    plot_data.rgb[2] = 1.0;

    GeneratePlotFromFunc(
        "plot.png", [](double x) { return x*x*x + x*x - 3*x - 3; }, 64, 1,
        3);
    reloadTexture(window, "plot.png");
    Vertex v{px, 0.f};

    v = plot_to_ortho_coords(v);

    points_->vertices.clear();
    points_->vertices.push_back(v);

    updateRenderObject(points_.get());
}


double f2(double x) {
    return sin(x) - (x-1)/2;
}

// Pochodna funkcji f'(x)
double df2(double x) {
    return cos(x) - 1/2;
}

double newton(double x0, double precision, int max_iter) {
    double x1;
    int iter = 0;
    do {
        if (df2(x0) == 0) {
            return x0;
        }
        x1 = x0 - f2(x0) / df2(x0);
        if (fabs(x1 - x0) < precision) break;
        x0 = x1;
        iter++;
    } while (iter < max_iter);

    return x1;
}

void on_key_s_pressed(GLFWwindow* window)
{
    double newton_start = 2.0;
    double precision = 1e-6;
    int max_iter = 100;

    plot_data.rgb[0] = 0.0;
    plot_data.rgb[1] = 1.0;
    plot_data.rgb[2] = 0.0;

    float px = (float)newton(newton_start, precision, max_iter);
    
    GeneratePlotFromFunc(
        "plot.png", [](double x) { return f2(x); }, 64, -10,
        10);
    reloadTexture(window, "plot.png");
    Vertex v{px, 0.f};

    v = plot_to_ortho_coords(v);

    points_->vertices.clear();
    points_->vertices.push_back(v);

    updateRenderObject(points_.get());
}

/////////////////////////////////////////////////////////////////////////
// Example of creating plot from click-points

void on_key_d_pressed(GLFWwindow* window)
{
    plot_data.plot_name = L"clicked";
    plot_data.line_type = L"dashed";
    plot_data.rgb[0] = 1.0;
    plot_data.rgb[1] = 0.0;
    plot_data.rgb[2] = 0.0;

    if (!GeneratePlotFromPoints(plot_filename_, plot_data.xs, plot_data.ys))
    {
        std::cerr << "Failed to generate initial plot image." << std::endl;
    }
    reloadTexture(window, "plot.png");

    // clear previous points
    points_->vertices.clear();
    plot_data.xs.clear();
    plot_data.ys.clear();

    updateRenderObject(points_.get());
}

/////////////////////////////////////////////////////////////////////////
// Empty slots for the students

void on_key_f_pressed(GLFWwindow* window) {
        const double PI = 3.141592653589793;

    plot_data.plot_name = L"sin/cos";
    plot_data.line_type = L"solid";
    plot_data.rgb[0] = 0.0;
    plot_data.rgb[1] = 0.0;
    plot_data.rgb[2] = 1.0;

    GenerateContinuousPlotFromFunc(
        "plot.png", [](double x) { return sinf((float)x); }, 64, -PI * 2.0,
        PI * 2.0);

    plot_data.rgb[0] = 1.0;
    plot_data.rgb[1] = 0.0;
    plot_data.rgb[2] = 0.0;

    GenerateContinuousPlotFromFunc(
        "plot.png", [](double x) { return cosf((float)x); }, 64, -PI * 2.0,
        PI * 2.0);

    reloadTexture(window, "plot.png");

    FinishContinuousPlot();

    points_->vertices.clear();
    plot_data.xs.clear();
    plot_data.ys.clear();
    updateRenderObject(points_.get());
}

bool solveGaussian(std::vector<std::vector<double>>& A, std::vector<double>& b, std::vector<double>& x)
{
    int n = A.size();

    for (int i = 0; i < n; i++)
    {
        // Wybór wiersza głównego
        double maxEl = fabs(A[i][i]);
        int maxRow = i;
        for (int k = i + 1; k < n; k++)
        {
            if (fabs(A[k][i]) > maxEl)
            {
                maxEl = fabs(A[k][i]);
                maxRow = k;
            }
        }

        // Zamiana
        std::swap(A[maxRow], A[i]);
        std::swap(b[maxRow], b[i]);

        // Eliminacja
        for (int k = i + 1; k < n; k++)
        {
            double c = -A[k][i] / A[i][i];
            for (int j = i; j < n; j++)
            {
                if (i == j) A[k][j] = 0;
                else A[k][j] += c * A[i][j];
            }
            b[k] += c * b[i];
        }
    }

    // Rozwiązanie układu
    x.assign(n, 0.0);
    for (int i = n - 1; i >= 0; i--)
    {
        x[i] = b[i] / A[i][i];
        for (int k = i - 1; k >= 0; k--)
        {
            b[k] -= A[k][i] * x[i];
        }
    }
    return true;
}


std::vector<double> polynomial_fit(const std::vector<double>& xs,
                                   const std::vector<double>& ys,
                                   int degree)
{
    int N = xs.size();
    int n = degree + 1;

    std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));
    std::vector<double> B(n, 0.0);

    std::vector<double> X(2 * degree + 1, 0.0);

    // obliczenie sum potęg x
    for (int i = 0; i < 2 * degree + 1; i++)
    {
        for (int k = 0; k < N; k++)
            X[i] += std::pow(xs[k], i);
    }

    // budowa macierzy A
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            A[i][j] = X[i + j];
    }

    // budowa wektora B
    for (int i = 0; i < n; i++)
    {
        for (int k = 0; k < N; k++)
            B[i] += std::pow(xs[k], i) * ys[k];
    }

    std::vector<double> coeffs;
    solveGaussian(A, B, coeffs);
    return coeffs;
}


// Oblicza y dla wielomianu
double poly_eval(const std::vector<double>& c, double x)
{
    double result = 0.0;
    double xn = 1.0;
    for (double ci : c)
    {
        result += ci * xn;
        xn *= x;
    }
    return result;
}

void on_key_g_pressed(GLFWwindow* window) {
        const auto& xs = plot_data.xs;
    const auto& ys = plot_data.ys;
    if (xs.empty() || ys.empty() || xs.size() != ys.size()) return;

    int maxDegree = 6;       
    const int samplesPerCurve = 200;

   
    points_->vertices.clear();

    
    float colors[][3] = {
        {1.0f, 0.0f, 0.0f}, 
        {0.0f, 1.0f, 0.0f}, 
        {0.0f, 0.0f, 1.0f}, 
        {1.0f, 1.0f, 0.0f}, 
        {1.0f, 0.0f, 1.0f}, 
        {0.0f, 1.0f, 1.0f}, 
    };

    
    double xmin = xs.front();
    double xmax = xs.back();
    if (xmin == xmax) { 
        xmin -= 1.0;
        xmax += 1.0;
    }
    double step = (xmax - xmin) / (samplesPerCurve - 1);

    
    for (int degree = 1; degree <= maxDegree; ++degree)
    {

        std::vector<double> coeffs = polynomial_fit(xs, ys, degree);

        int ci = (degree - 1) % (sizeof(colors)/sizeof(colors[0]));
        plot_data.rgb[0] = colors[ci][0];
        plot_data.rgb[1] = colors[ci][1];
        plot_data.rgb[2] = colors[ci][2];

        for (int i = 0; i < samplesPerCurve; ++i)
        {
            double x = xmin + i * step;
            double y = poly_eval(coeffs, x);

            Vertex v;
            v.x = (float)x;
            v.y = (float)y;
            v = plot_to_ortho_coords(v);

            points_->vertices.push_back(v);
        }


        updateRenderObject(points_.get());

    }

    plot_data.rgb[0] = 1.0;
    plot_data.rgb[1] = 1.0;
    plot_data.rgb[2] = 1.0;
    for (size_t i = 0; i < xs.size(); ++i)
    {
        Vertex pv;
        pv.x = (float)xs[i];
        pv.y = (float)ys[i];
        pv = plot_to_ortho_coords(pv);
        points_->vertices.push_back(pv);
    }
    updateRenderObject(points_.get());
}


void fft(CArray& x) {
    const size_t N = x.size();
    if (N <= 1) return;

    // Podział na parzyste i nieparzyste
    CArray even = x[std::slice(0, N / 2, 2)];
    CArray odd = x[std::slice(1, N / 2, 2)];

    fft(even);
    fft(odd);

    for (size_t k = 0; k < N / 2; ++k) {
        Complex t = std::polar(1.0, -2 * PI * k / N) * odd[k];
        x[k] = even[k] + t;
        x[k + N / 2] = even[k] - t;
    }
}

void on_key_h_pressed(GLFWwindow* window) {
    const int N = 1024;           
    const double Fs = 1000.0;      
    
    CArray data(N);
    for (int i = 0; i < N; ++i) {
        double t = (double)i / Fs;
        // val: funkcja sygnału wejściowego
        // Sygnał: suma 50Hz (amp 1.0) i 120Hz (amp 0.5)
        double val = 1.0 * sin(2 * PI * 50 * t) + 0.5 * sin(2 * PI * 120 * t);
        data[i] = Complex(val, 0.0);
    }

    fft(data);

    std::vector<double> freqs;
    std::vector<double> amplitudes;
    
    for (int i = 0; i < N / 2; ++i) {
        double freq = i * Fs / N; 
        double amp = std::abs(data[i]) / (N / 2); 
        
        freqs.push_back(freq);
        amplitudes.push_back(amp);
    }

    plot_data.plot_name = L"Widmo FFT (Amplituda)";
    plot_data.range_x_min = 0.0f;
    plot_data.range_x_max = (float)(Fs / 2.0); 
    plot_data.range_y_min = 0.0f;
    plot_data.range_y_max = 1.2f; 
    
    plot_data.rgb[0] = 0.0f; plot_data.rgb[1] = 1.0f; plot_data.rgb[2] = 0.0f;


    GeneratePlotFromPoints("plot.png", freqs, amplitudes);
    reloadTexture(window, "plot.png");

    points_->vertices.clear();
    updateRenderObject(points_.get());
}
void on_key_m_pressed(GLFWwindow* window) {
    std::vector<fs::path> mp3_files;
    std::cout << "\n--- Wybierz plik MP3 z listy ---\n";
    int idx = 0;
    for (const auto& entry : fs::directory_iterator(".")) {
        if (entry.path().extension() == ".mp3") {
            mp3_files.push_back(entry.path());
            std::cout << idx++ << ": " << entry.path().filename() << "\n";
        }
    }

    if (mp3_files.empty()) {
        std::cerr << "Nie znaleziono plików MP3 w bieżącym folderze!\n";
        return;
    }

    int choice;
    std::cout << "Podaj numer pliku: ";
    std::cin >> choice;
    if (choice < 0 || choice >= mp3_files.size()) return;

    drmp3 mp3;
    if (!drmp3_init_file(&mp3, mp3_files[choice].string().c_str(), NULL)) {
        std::cerr << "Błąd podczas otwierania pliku!\n";
        return;
    }

    // Pobieramy fragment sygnału. Musi być długości 2^N 
    const size_t N = 16384; 
    std::vector<float> pSampleData(N * mp3.channels);
    drmp3_uint64 framesRead = drmp3_read_pcm_frames_f32(&mp3, N, pSampleData.data());
    
    double Fs = mp3.sampleRate;
    drmp3_uninit(&mp3);


    CArray complexData(N);
    for (size_t i = 0; i < N; ++i) {
        complexData[i] = Complex(pSampleData[i * mp3.channels], 0.0);
    }


    fft(complexData);

    std::vector<double> freqs;
    std::vector<double> amplitudes;
    double max_amp = 0;

    for (size_t i = 0; i < N / 2; ++i) {
        double freq = i * Fs / N;
        double amp = std::abs(complexData[i]);
        if (amp > max_amp) max_amp = amp;
        
        freqs.push_back(freq);
        amplitudes.push_back(amp);
    }

    plot_data.plot_name = L"Widmo MP3: " + mp3_files[choice].filename().wstring();
    plot_data.range_x_min = 0.0f;
    plot_data.range_x_max = 5000.0f; 
    plot_data.range_y_min = 0.0f;
    plot_data.range_y_max = (float)max_amp * 1.1f;
    
    plot_data.rgb[0] = 1.0f; plot_data.rgb[1] = 0.5f; plot_data.rgb[2] = 0.0f;

    GeneratePlotFromPoints("plot.png", freqs, amplitudes);
    reloadTexture(window, "plot.png");
    
    std::cout << "Zakończono analizę widmową pliku.\n";
}

// clear user interaction buffers
void on_key_clear(GLFWwindow* window)
{
    // clear previous points
    points_->vertices.clear();
    plot_data.xs.clear();
    plot_data.ys.clear();
    updateRenderObject(points_.get());

    lines_x_->vertices.clear();
    plot_data.user_range_x[0] = 0.0;
    plot_data.user_range_x[1] = 0.0;
    updateRenderObject(lines_x_.get());

    lines_y_->vertices.clear();
    plot_data.user_range_y[0] = 0.0;
    plot_data.user_range_y[1] = 0.0;
    updateRenderObject(lines_y_.get());
}

/////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
    // --- initialize from command line ---
    plot_data.pad_x = plot_data.pix_x / 20;
    plot_data.pad_y = plot_data.pix_y / 20;

    // --- parse command line options ---
    CLI::App app{ "Numerical Recipes Plotter" };

    app.add_option("-x, --plot-width", plot_data.pix_x, "Plot width in pixels")
        ->default_val(640);
    app.add_option("-y, --plot-height", plot_data.pix_y,
                   "Plot height in pixels")
        ->default_val(480);

    app.add_option("-a,--pad-width", plot_data.pad_x,
                   "Plot padding width in pixels")
        ->default_val(32);
    app.add_option("-s,--pad-height", plot_data.pad_y,
                   "Plot padding height in pixels")
        ->default_val(24);

    app.add_option("-z,--range-minx", plot_data.range_x_min, "Plot min X")
        ->default_val(-10.f);
    app.add_option("-c,--range-maxx", plot_data.range_x_max, "Plot max X")
        ->default_val(10.f);

    app.add_option("-t,--range-miny", plot_data.range_y_min, "Plot min Y")
        ->default_val(-10.f);
    app.add_option("-u,--range-maxy", plot_data.range_y_max, "Plot max Y")
        ->default_val(10.f);

    // 2. Parse the command line.
    // This macro includes a try/catch block and will exit cleanly on --help.
    CLI11_PARSE(app, argc, argv);


    if (!GenerateEmptyPlot(plot_filename_))
    {
        std::cerr << "Failed to generate initial plot image." << std::endl;
        return -1;
    }

    // --- GLFW initialization ---
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    // Load initial image dimensions to create a well-sized window
    int initial_width, initial_height, nrChannels;
    unsigned char* initial_data =
        stbi_load(plot_filename_.c_str(), &initial_width, &initial_height,
                  &nrChannels, 0);
    if (!initial_data)
    {
        std::cerr << "Failed to load initial texture to determine size."
                  << std::endl;
        glfwTerminate();
        return -1;
    }
    stbi_image_free(initial_data);

    // --- GLFW window creation ---
    GLFWwindow* window =
        glfwCreateWindow(initial_width, initial_height,
                         "PNG Viewer | Press 'R' to refresh", NULL, NULL);
    if (window == NULL)
    {
        std::cerr << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);

    // --- Set Callbacks ---
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetKeyCallback(window, key_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetCursorPosCallback(window, cursor_position_callback);

    // --- GLAD ---
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cerr << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    {
        float pcol[] = { 1.f, 0.f, 0.f };
        points_.reset(createRenderObject(GL_POINTS, pcol, 4.f));
    }

    {
        float lcol[] = { 0.f, 1.f, 0.f };
        lines_x_.reset(createRenderObject(GL_LINES, lcol, 2.f));
    }

    {
        float lcol[] = { 0.f, 0.f, 1.f };
        lines_y_.reset(createRenderObject(GL_LINES, lcol, 2.f));
    }

    // --- Shaders (compile and link) ---
    unsigned int shaderProgram = 0, overlayProgram = 0;
    {
        unsigned int vertexShader = glCreateShader(GL_VERTEX_SHADER);
        glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
        glCompileShader(vertexShader);
        unsigned int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
        glCompileShader(fragmentShader);
        shaderProgram = glCreateProgram();
        glAttachShader(shaderProgram, vertexShader);
        glAttachShader(shaderProgram, fragmentShader);
        glLinkProgram(shaderProgram);
        glDeleteShader(vertexShader);
        glDeleteShader(fragmentShader);
    }

    {
        unsigned int vertexShader = glCreateShader(GL_VERTEX_SHADER);
        glShaderSource(vertexShader, 1, &vertexOvlSource, NULL);
        glCompileShader(vertexShader);
        unsigned int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(fragmentShader, 1, &fragmentOvlSource, NULL);
        glCompileShader(fragmentShader);
        overlayProgram = glCreateProgram();
        glAttachShader(overlayProgram, vertexShader);
        glAttachShader(overlayProgram, fragmentShader);
        glLinkProgram(overlayProgram);
        glDeleteShader(vertexShader);
        glDeleteShader(fragmentShader);
    }


    // --- Vertex Data for a Plot Quad ---
    float vertices[] = { // positions      // texture coords
                         1.0f, 1.0f,  0.0f, 1.0f,  0.0f,  1.0f, -1.0f,
                         0.0f, 1.0f,  1.0f, -1.0f, -1.0f, 0.0f, 0.0f,
                         1.0f, -1.0f, 1.0f, 0.0f,  0.0f,  0.0f
    };
    unsigned int indices[] = { 0, 1, 3, 1, 2, 3 };
    unsigned int VBO, VAO, EBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices,
                 GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float),
                          (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float),
                          (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    // --- Texture Generation ---
    glGenTextures(1, &texture_id_); // Use the global texture ID
    glBindTexture(GL_TEXTURE_2D, texture_id_);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    // Load the initial texture data using our new function
    reloadTexture(window, plot_filename_);

    // --- Render loop ---
    while (!glfwWindowShouldClose(window))
    {
        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        // --- Dynamic Viewport for Aspect Ratio Correction ---
        int win_width, win_height;
        glfwGetFramebufferSize(window, &win_width, &win_height);
        float texture_aspect = (texture_width_ > 0 && texture_height_ > 0)
            ? (float)texture_width_ / (float)texture_height_
            : 1.0f;
        int view_width = win_width;
        int view_height = (int)(view_width / texture_aspect);
        if (view_height > win_height)
        {
            view_height = win_height;
            view_width = (int)(view_height * texture_aspect);
        }
        int view_x = (win_width - view_width) / 2;
        int view_y = (win_height - view_height) / 2;
        glViewport(view_x, view_y, view_width, view_height);

        // --- Draw Quad ---
        glUseProgram(shaderProgram);
        float projection[4][4] = { { 1.0f, 0.0f, 0.0f, 0.0f },
                                   { 0.0f, 1.0f, 0.0f, 0.0f },
                                   { 0.0f, 0.0f, 1.0f, 0.0f },
                                   { 0.0f, 0.0f, 0.0f, 1.0f } };
        glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "projection"), 1,
                           GL_FALSE, &projection[0][0]);

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, texture_id_);
        glBindVertexArray(VAO);
        glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

        // --- Draw Overlays ---
        glUseProgram(overlayProgram);
        glUniformMatrix4fv(glGetUniformLocation(overlayProgram, "projection"),
                           1, GL_FALSE, &projection[0][0]);
        drawRenderObject(points_.get(), overlayProgram);
        drawRenderObject(lines_x_.get(), overlayProgram);
        drawRenderObject(lines_y_.get(), overlayProgram);

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // --- Cleanup ---
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteBuffers(1, &EBO);
    glDeleteProgram(shaderProgram);
    glDeleteTextures(1, &texture_id_);

    glfwTerminate();
    return 0;
}

/////////////////////////////////////////////////////////////////////////
/**
 * @brief Reloads the PNG from disk and updates the active OpenGL texture.
 * @param filename The path to the PNG image file.
 */
void reloadTexture(GLFWwindow* window, const std::string& filename)
{
    int new_width, new_height, nrChannels;
    // Force loading the image with 4 channels (RGBA) for consistency
    unsigned char* data =
        stbi_load(filename.c_str(), &new_width, &new_height, &nrChannels, 4);

    if (data)
    {
        texture_width_ = new_width;
        texture_height_ = new_height;

        glBindTexture(GL_TEXTURE_2D, texture_id_);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, texture_width_, texture_height_,
                     0, GL_RGBA, GL_UNSIGNED_BYTE, data);
        glGenerateMipmap(GL_TEXTURE_2D);

        glfwSetWindowSize(window, texture_width_, texture_height_);

        stbi_image_free(data);
        std::cout << "Texture '" << filename << "' reloaded (" << texture_width_
                  << "x" << texture_height_ << ")." << std::endl;
    }
    else
    {
        std::cerr << "Failed to reload texture from '" << filename << "'."
                  << std::endl;
    }
}

/////////////////////////////////////////////////////////////////////////
/**
 * @brief Handles all key press events for the application.
 */
void key_callback(GLFWwindow* window, int key, int scancode, int action,
                  int mods)
{
    // We only care about key press events, not releases
    if (action != GLFW_PRESS)
    {
        return;
    }

    key_pressed_ = key;

    switch (key)
    {
        case GLFW_KEY_ESCAPE: glfwSetWindowShouldClose(window, true); break;
        case GLFW_KEY_R: reloadTexture(window, "plot.png"); break;
        case GLFW_KEY_A: on_key_a_pressed(window); break;
        case GLFW_KEY_S: on_key_s_pressed(window); break;
        case GLFW_KEY_D: on_key_d_pressed(window); break;
        case GLFW_KEY_F: on_key_f_pressed(window); break;
        case GLFW_KEY_G: on_key_g_pressed(window); break;
        case GLFW_KEY_Z: on_key_z_pressed(window); break;
        case GLFW_KEY_H: on_key_h_pressed(window); break;
        case GLFW_KEY_M: on_key_m_pressed(window); break;

        case GLFW_KEY_C: on_key_clear(window); break;
    }
}

/////////////////////////////////////////////////////////////////////////

void mouse_button_x_range(const double& xpos, const double& ypos)
{
    if (lines_x_->vertices.size() <= 2)
    {
        Vertex v0 = { (float)(xpos / plot_data.pix_x) * 2.f - 1.f, -1.f };
        Vertex v1 = { (float)(xpos / plot_data.pix_x) * 2.f - 1.f, 1.f };

        lines_x_->vertices.push_back(v0);
        lines_x_->vertices.push_back(v1);

        float plot_x = float(xpos - plot_data.pad_x)
            / (plot_data.pix_x - plot_data.pad_x * 2);
        plot_x = plot_data.range_x_min
            + plot_x * (plot_data.range_x_max - plot_data.range_x_min);
        plot_data.user_range_x[0] = plot_x;
    }
    else
    {
        float plot_x = float(xpos - plot_data.pad_x)
            / (plot_data.pix_x - plot_data.pad_x * 2);
        plot_x = plot_data.range_x_min
            + plot_x * (plot_data.range_x_max - plot_data.range_x_min);
        plot_data.user_range_x[1] = plot_x;

        if (plot_data.user_range_x[1] < plot_data.user_range_x[0])
            std::swap(plot_data.user_range_x[1], plot_data.user_range_x[0]);

        key_pressed_ = 0;
    }
    updateRenderObject(lines_x_.get());
}

/////////////////////////////////////////////////////////////////////////

void mouse_button_y_range(const double& xpos, const double& ypos)
{
    if (lines_y_->vertices.size() <= 2)
    {
        Vertex v0 = { -1.f, float(1.f - ypos / plot_data.pix_y) * 2.f - 1.f };
        Vertex v1 = { 1.f, float(1.f - ypos / plot_data.pix_y) * 2.f - 1.f };

        lines_y_->vertices.push_back(v0);
        lines_y_->vertices.push_back(v1);

        float plot_y = 1.f
            - (float(ypos - plot_data.pad_y)
               / (plot_data.pix_y - plot_data.pad_y * 2));
        plot_y = plot_data.range_y_min
            + plot_y * (plot_data.range_y_max - plot_data.range_y_min);
        plot_data.user_range_y[0] = plot_y;
    }
    else
    {
        float plot_y = 1.f
            - (float(ypos - plot_data.pad_y)
               / (plot_data.pix_y - plot_data.pad_y * 2));
        plot_y = plot_data.range_y_min
            + plot_y * (plot_data.range_y_max - plot_data.range_y_min);
        plot_data.user_range_y[1] = plot_y;

        if (plot_data.user_range_y[1] < plot_data.user_range_y[0])
            std::swap(plot_data.user_range_y[1], plot_data.user_range_y[0]);

        key_pressed_ = 0;
    }
    updateRenderObject(lines_y_.get());
}

/////////////////////////////////////////////////////////////////////////
/**
 * @brief Handles mouse button events, storing left-click coordinates in a
 * buffer.
 */
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
    {
        double xpos, ypos;
        glfwGetCursorPos(window, &xpos, &ypos);

        if (key_pressed_ == GLFW_KEY_X)
        {
            mouse_button_x_range(xpos, ypos);
        }
        else if (key_pressed_ == GLFW_KEY_Y)
        {
            mouse_button_y_range(xpos, ypos);
        }
        else
        {
            // update plot resources
            {
                Vertex v = pixel_to_plot({(float)xpos, (float)ypos});
                plot_data.xs.push_back(v.x);
                plot_data.ys.push_back(v.y);
            }

            // update rendering resources
            {
                Vertex v = pixel_to_ortho_coords({(float)xpos, (float)ypos});

                // correct y due to inverted pixel origin from glfw
                v.y=-v.y;
                points_->vertices.push_back(v);

                updateRenderObject(points_.get());
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////////

void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
{
    if (key_pressed_ == GLFW_KEY_X)
    {
        Vertex v0 = { float(xpos / plot_data.pix_x) * 2.f - 1.f, -1.f };
        Vertex v1 = { float(xpos / plot_data.pix_x) * 2.f - 1.f, 1.f };
        uint8_t i0 = 0, i1 = 1;

        if (lines_x_->vertices.size() <= 2)
        {
            lines_x_->vertices.resize(2);
        }
        else
        {
            i0 = 2;
            i1 = 3;
        }
        lines_x_->vertices[i0] = v0;
        lines_x_->vertices[i1] = v1;
        updateRenderObject(lines_x_.get());
    }
    else if (key_pressed_ == GLFW_KEY_Y)
    {
        Vertex v0 = { -1.f, float(1.f - ypos / plot_data.pix_y) * 2.f - 1.f };
        Vertex v1 = { 1.f, float(1.f - ypos / plot_data.pix_y) * 2.f - 1.f };
        uint8_t i0 = 0, i1 = 1;

        if (lines_y_->vertices.size() <= 2)
        {
            lines_y_->vertices.resize(2);
        }
        else
        {
            i0 = 2;
            i1 = 3;
        }
        lines_y_->vertices[i0] = v0;
        lines_y_->vertices[i1] = v1;
        updateRenderObject(lines_y_.get());
    }
}

/////////////////////////////////////////////////////////////////////////
/**
 * @brief Callback for window resize events. Kept for completeness.
 */
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    // The viewport is handled dynamically in the main render loop to maintain
    // aspect ratio. This function is still required by GLFW but can be left
    // empty for our purpose.
}