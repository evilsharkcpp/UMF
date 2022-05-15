#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>

using namespace std;

//Lab4

class BCG
{
    vector<vector<double>> a;
    vector<double> b;
    vector<double> x;
    vector<double> r;
    vector<double> rv;
    vector<double> v;
    vector<double> p;
    int maxIterations;

    double norm(vector<double>& x)
    {
        double res = 0;
        for (int i = 0; i < x.size(); i++)
            res += x[i] * x[i];
        return sqrt(res);
    }
    double normv(vector<double>& x)
    {

    }

    void multMatrixByVector(vector<double>& res, vector<vector<double>>& mat, vector<double>& v)
    {
        for (int i = 0; i < res.size(); i++)
        {
            double sum = 0;
            for (int j = 0; j < mat[i].size(); j++)
                sum += mat[i][j] * v[j];
            res[i] = sum;
        }
    }

    double scalar(vector<double>& a, vector<double>& b)
    {
        double sum = 0;
        for (int i = 0; i < a.size(); i++)
            sum += a[i] * b[i];
        return sum;
    }

    vector<vector<double>> aT()
    {
        vector<vector<double>> res(a[0].size(), vector<double>(a.size()));
        for (int i = 0; i < a.size(); i++)
            for (int j = 0; j < a[i].size(); j++)
                res[j][i] = a[i][j];
        return res;
    }
public:

    BCG(vector<vector<double>>& mat, vector<double>& f, int iterations)
    {
        a = vector<vector<double>>(mat);
        b = vector<double>(f);
        x.resize(b.size());
        rv.resize(b.size());
        r.resize(b.size());
        v.resize(b.size());
        p.resize(b.size());
        maxIterations = iterations;
    }

    void calculate(vector<double>& x_start, double eps = 1e-8)
    {
       //Init
        vector<double> xk(x_start);
        multMatrixByVector(r, a, x_start);
        for (int i = 0; i < r.size(); i++)
        {
            r[i] = b[i] - r[i];
            rv[i] = r[i];
        }
        double ro = 1, alpha = 1, omega = 1;
        for (int i = 0; i < p.size(); i++)
        {
            p[i] = 0;
            v[i] = 0;
        }
        //calc
        double disp = 1;
        int i = 0;
        while (i < maxIterations && disp > eps)
        {
            double rok = scalar(rv, r);
            double beta = (rok / ro) * (alpha / omega);
            for (int i = 0; i < p.size(); i++)
                p[i] = r[i] + beta * (p[i] - omega * v[i]);
            multMatrixByVector(v, a, p);
            alpha = rok / scalar(rv, v);
            vector<double> s(p.size());
            for (int i = 0; i < s.size(); i++)
                s[i] = r[i] - alpha * v[i];
            vector<double> t(p.size());
            multMatrixByVector(t, a, s);
            omega = scalar(t, s) / scalar(t, t);
            for (int i = 0; i < x.size(); i++)
            {
                x[i] = xk[i] + omega * s[i] + alpha * p[i];
                xk[i] = x[i];
                r[i] = s[i] - omega * t[i];
            }
            ro = rok;
            vector<double> f(b.size());
            multMatrixByVector(f, a, x);
            for (int i = 0; i < f.size(); i++)
                f[i] -= b[i];
            disp = norm(f);
        }
    }

    void getResult(vector<double>& res)
    {
        for (int i = 0; i < x.size(); i++)
            res[i] = x[i];
    }
};

//Lab3

struct Element
{
    int number;
    vector<int> nodes;
    Element()
    {
        number = 0;
        nodes.resize(4);
    }
    void setNodes(int i, int kj)
    {
        nodes[0] = number + i;
        nodes[1] = number + i + 1;
        nodes[2] = number + kj + i + 1;
        nodes[3] = number + kj + i + 2;
    }
    int getNode(int num)
    {
        return nodes[num];
    }
};

struct Point
{
    double X;
    double Y;
    Point()
    {
        X = 0;
        Y = 0;
    }
    Point(double x, double y)
    {
        X = x;
        Y = y;
    }
};

class Matrix
{
    //vector<vector<double>> a;
public:
    vector<vector<double>> a;
    Matrix()
    {

    }
    Matrix(int n)
    {
        a.resize(n, vector<double>(n));
    }
    double getElement(int i, int j)
    {
        return a[i][j];
    }
    void setElement(int i, int j, double value)
    {
        a[i][j] += value;
    }
    void print()
    {
        ofstream out("out.txt");
        for (int i = 0; i < a.size(); i++)
        {
            for (int j = 0; j < a[i].size(); j++)
            {
                out << setprecision(5) << fixed << a[i][j] << " ";
            }
            out << endl;
        }
    }
};

class FEM
{
    vector<vector<Point>> grid;
    vector<Element> elements;
    Matrix a;
    vector<double> b;
    vector<double> q;

    double sigma;
    double hi;
    double lambda;
    double omega;
    vector<vector<int>> M = {
        {4,2,2,1},
        {2,4,1,2},
        {2,1,4,2},
        {1,2,2,4}
    };
    vector<vector<int>> G = {
        {2,-2,1,-1},
        {-2,2,-1,1},
        {1,-1,2,-2},
        {-1,1,-2,2}
    };

    template<typename T, typename G>
    void multMatrixByConst(vector<vector<T>>& res, vector<vector<G>>& mat, double c)
    {
        for (int i = 0; i < res.size(); i++)
        {
            for (int j = 0; j < res[i].size(); j++)
            {
                res[i][j] = c * mat[i][j];
            }
        }
    }

    template<typename T, typename G, typename Q>
    void AddMatrixs(vector<vector<T>>& res, vector<vector<G>>& a, vector<vector<Q>>& b)
    {
        for (int i = 0; i < res.size(); i++)
            for (int j = 0; j < res[i].size(); j++)
                res[i][j] = a[i][j] + b[i][j];
    }
    template<typename T, typename G, typename Q>
    void SubMatrixs(vector<vector<T>>& res, vector<vector<G>>& a, vector<vector<Q>>& b)
    {
        for (int i = 0; i < res.size(); i++)
            for (int j = 0; j < res[i].size(); j++)
                res[i][j] = a[i][j] - b[i][j];
    }

    int getLocalNum(int j, int nx, int num)
    {
        int i = num / nx;
        return num - i * nx; // return local i by num, total elem in grid and local j
    }

    void addLocalInGlobal(int gi, int gj, vector<vector<double>>& local, vector<double>& blocal)
    {
        int shiftI = (gi == 0) ? 0 : 1;
        int shiftJ = (gj == 0) ? 0 : 1;
        a.setElement(2 * gi, 2 * gj, local[0][0]);
        a.setElement(2 * gi, 2 * gj + 1, local[0][1]);
        a.setElement(2 * gi + 1, 2 * gj, local[1][0]);
        a.setElement(2 * gi + 1, 2 * gj + 1, local[1][1]);
        b[2 * gi] += blocal[0];
        b[2 * gi + 1] += blocal[1];
    }

   
public:

    double fs(double x, double y)
    {
        return 0;
    }

    double fc(double x, double y)
    {
        return 0;
    }

    void Init(int nx, int ny, double ax = 0, double bx = 1, double ay = 0, double by = 1)
    {
        elements.resize(nx * ny);
        double hx = (bx - ax) / nx;
        double hy = (by - ay) / ny;

        for (int i = 0; i <= ny; i++)
        {
            vector<Point> tmp;
            for (int j = 0; j < nx; j++)
            {
                tmp.push_back(Point(ax + hx * j, ay + hy * i));

            }
            tmp.push_back(Point(bx, ay + hy * i));
            grid.push_back(tmp);
        }
        int num = 0;
        for (int i = 0; i < grid.size() - 1; i++)
            for (int j = 0; j < grid[i].size() - 1; j++)
            {
                elements[num].number = num;
                elements[num++].setNodes(i, nx);
            }
        a = Matrix(2 * grid.size() * grid[0].size());
        b.resize(2 * grid.size() * grid[0].size());
        q.resize(2 * grid.size() * grid[0].size());
    }

    void setParams(double sig, double omeg, double khi, double lam)
    {
        sigma = sig;
        omega = omeg;
        hi = khi;
        lambda = lam;
    }

    void setMatrix()
    {
        vector<vector<double>> aij(2, vector<double>(2));
        vector<double> bij(2);
        vector<vector<double>> pij(M.size(), vector<double>(M.size()));
        vector<vector<double>> tmp(M.size(), vector<double>(M.size()));
        vector<vector<double>> cij(M.size(), vector<double>(M.size()));
        vector<vector<double>> vec(M.size(), vector<double>(M.size()));
        int nx = grid[0].size();
        for (int k = 0; k < elements.size(); k++)
        {
            double hx = abs(grid[elements[k].getNode(0) / nx][getLocalNum(0, nx, elements[k].getNode(0))].X - grid[elements[k].getNode(1) / nx][getLocalNum(1, nx, elements[k].getNode(1))].X);
            double hy = abs(grid[elements[k].getNode(0) / nx][getLocalNum(0, nx, elements[k].getNode(0))].Y - grid[elements[k].getNode(2) / nx][getLocalNum(2, nx, elements[k].getNode(2))].Y);
            double coef = omega * sigma * hx * hy / 36.0;

            multMatrixByConst(cij, M, coef);
            multMatrixByConst(vec, M, hx * hy / 36.0);
            multMatrixByConst(pij, G, lambda * hx / (6 * hy));
            multMatrixByConst(tmp, G, lambda * hy / (6 * hx));
            AddMatrixs(pij, pij, tmp);
            multMatrixByConst(tmp, M, hi * omega * omega);
            SubMatrixs(pij, pij, tmp);
            for (int i = 0; i < 4; i++)
            {
                int gi = elements[k].getNode(i);
                for (int j = 0; j < 4; j++)
                {
                    int gj = elements[k].getNode(j);
                    aij[0][0] = pij[i][j];
                    aij[0][1] = -cij[i][j];
                    aij[1][0] = cij[i][j];
                    aij[1][1] = pij[i][j];
                    int lj = getLocalNum(j, nx, gj);
                    int li = gj / nx;
                    double x = grid[li][lj].X;
                    double y = grid[li][lj].Y;
                    bij[0] = fs(x, y) * vec[i][j];
                    bij[1] = fc(x, y) * vec[i][j];
                    addLocalInGlobal(gi, gj, aij, bij);
                }
            }
        }
        a.print();
    }

    double GetBound(double x, double y, int num, int type)
    {
        if (!type)
        {
            switch (num)
            {
            case 1:
                return x;
                break;
            case 2:
                return 0;
                break;
            case 3:
                return x;
                break;
            case 4:
                return 1;
                break;
            }
        }
        else
        {

        }
    }

    void FirstBound()
    {
        ifstream in("bound.txt");
        int num;
        while (in >> num)
        {
            switch (num)
            {
            case 1:
                for (int i = 0; i < grid[0].size(); i++)
                {
                    for (int j = 0; j < a.a[i].size(); j++)
                    {
                        a.a[2 * i][j] = 0;
                        a.a[2 * i + 1][j] = 0;
                    }
                    a.a[2 * i][2 * i] = 1;
                    a.a[2 * i + 1][2 * i + 1] = 1;
                    b[2 * i] = GetBound(grid[0][i].X, grid[0][i].Y, num, 0);
                    b[2 * i + 1] = GetBound(grid[0][i].X, grid[0][i].Y, num, 1);

                }
                break;
            case 2:
                for (int i = 0; i < grid.size()*grid[0].size(); i += grid[0].size())
                {
                    for (int j = 0; j < a.a[i].size(); j++)
                    {
                        a.a[2 * i][j] = 0;
                        a.a[2 * i + 1][j] = 0;
                    }
                    a.a[2 * i][2 * i] = 1;
                    a.a[2 * i + 1][2 * i + 1] = 1;
                    int q = i / grid[0].size();
                    b[2 * i] = GetBound(grid[q][0].X, grid[q][0].Y, num, 0);
                    b[2 * i + 1] = GetBound(grid[q][0].X, grid[q][0].Y, num, 1);
                }
                break;
            case 3:
                for (int k = 0; k < grid[0].size(); k++)
                {
                    int i = (grid.size() - 1) * grid[0].size() + k;
                    for (int j = 0; j < a.a[i].size(); j++)
                    {
                        a.a[2 * i][j] = 0;
                        a.a[2 * i + 1][j] = 0;
                    }
                    a.a[2 * i][2 * i] = 1;
                    a.a[2 * i + 1][2 * i + 1] = 1;
                    b[2 * i] = GetBound(grid[grid.size() - 1][k].X, grid[grid.size() - 1][k].Y, num, 0);
                    b[2 * i + 1] = GetBound(grid[grid.size() - 1][k].X, grid[grid.size() - 1][k].Y, num, 1);
                }
                break;
            case 4:
                for (int i = grid[0].size() - 1; i < grid.size() * grid[0].size() + grid[0].size() - 1; i += grid[0].size())
                {
                    for (int j = 0; j < a.a[i].size(); j++)
                    {
                        a.a[2 * i][j] = 0;
                        a.a[2 * i + 1][j] = 0;
                    }
                    a.a[2 * i][2 * i] = 1;
                    a.a[2 * i + 1][2 * i + 1] = 1;
                    int q = i / grid[0].size();
                    b[2 * i] = GetBound(grid[q][grid[0].size() - 1].X, grid[q][grid[0].size() - 1].Y, num, 0);
                    b[2 * i + 1] = GetBound(grid[q][grid[0].size() - 1].X, grid[q][grid[0].size() - 1].Y, num, 1);
                }
                break;
            }
        }
        a.print();
    }
    
    void solve()
    {
        BCG solv(a.a, b, 2);
        auto x = vector<double>(b.size(), 0);
        solv.calculate(x);
        solv.getResult(q);
    }
    
    void print()
    {
        for (int i = 0; i < grid.size(); i++)
            for (int j = 0; j < grid[i].size(); j++)
            {
                int k = i * grid[i].size() + j;
                cout << q[2 * k] << endl;
            }
            
    }
};

int main()
{
    vector<vector<double>> a = {
        {1,2,3},
        {4,5,6},
        {7,8,9}
    };
    vector<double> b = {
        {14,32,50}
    };
    BCG bcg(a, b, 100);
    auto x = vector<double>(3, 0);
    //bcg.calculate(x);
    FEM fem;
    fem.Init(2, 2);
    fem.setParams(0, 0, 0, 1);
    fem.setMatrix();
    fem.FirstBound();
    fem.solve();
    fem.print();
}