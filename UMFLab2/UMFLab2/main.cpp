#include <iostream>
#include <vector>
#include <fstream>
#include <string>

using namespace std;

class LUMethod
{
    vector<double> al;
    vector<double> au;
    vector<double> di;
    vector<double> b;
    vector<double> q;
    vector<double> y;

    double getElement(int i, int j)
    {
        if (abs(j - i) == 1 || abs(j - i) == 0)
        {
            switch (j - i)
            {
            case -1:
                return al[j + 1];
                break;
            case 0:
                return di[j];
                break;
            case 1:
                return au[j - 1];
                break;
            }
        }
        return 0;
    }

    void setElement(int i, int j, double elem)
    {
        if (abs(j - i) == 1 || abs(j - i) == 0)
        {
            switch (j - i)
            {
            case -1:
                al[j + 1] = elem;
                break;
            case 0:
                di[j] = elem;
                break;
            case 1:
                au[j - 1] = elem;
                break;
            }
        }
    }
public:
    LUMethod(vector<vector<double>>& a, vector<double>& b)
    {
        if (a.size() == 3)
        {
            al.resize(a[0].size());
            au.resize(a[0].size());
            di.resize(a[0].size());
            this->b.resize(a[0].size());
            q.resize(a[0].size());
            for (int i{ 0 }; i < a[0].size(); i++)
            {
                al[i] = a[0][i];
                di[i] = a[1][i];
                au[i] = a[2][i];
                this->b[i] = b[i];
                y.resize(b.size());
            }
        }
        else throw string("wrong matrix");
    }

    void LU()
    {
        for (int i{ 0 }; i < di.size(); i++)
        {
            for (int j{ 0 }; j < di.size(); j++)
            {
                if (j >= i)
                {
                    double sum{ 0 };
                    for (int k{ 0 }; k < i; k++)
                        sum += getElement(i, k) * getElement(k, j);
                    setElement(i, j, getElement(i, j) - sum);

                }
                else
                {
                    double sum{ 0 };
                    for (int k{ 0 }; k < j; k++)
                        sum += getElement(i, k) * getElement(k, j);
                    setElement(i, j, (getElement(i, j) - sum) / getElement(j, j));
                }
            }
        }
    }

    void forwardStep()
    {
        for (int i = 0; i < di.size(); i++)
        {
            double sum = 0;
            for (int j = 0; j < i; j++)
                sum += getElement(i, j) * y[j];
            y[i] = (b[i] - sum);
        }
    }

    void backStep()
    {
        for (int i = di.size() - 1; i >= 0; i--)
        {
            double sum = 0;
            for (int j = i + 1; j < di.size(); j++)
                sum += getElement(i, j) * q[j];
            q[i] = (y[i] - sum) / getElement(i, i);
        }
    }
    
    void getQ(vector<double>& a)
    {
        for (int i = 0; i < q.size(); i++)
            a[i] = q[i];
    }
};

struct Element
{
    double lambda;
    double gamma;
};

class Matrix
{

public:
    vector<vector<double>> a;

    Matrix()
    {

    }

    Matrix(int n)
    {
        a.resize(3);
        for (size_t i{ 0 }; i < 3; i++)
            a[i].resize(n);
    }

    double getElement(int i, int j)
    {
        if (abs(j - i) == 1 || abs(j - i) == 0)
        {
            switch (j - i)
            {
            case -1:
                return a[0][j + 1];
                break;
            case 0:
                return a[1][j];
                break;
            case 1:
                return a[2][j - 1];
                break;
            }
        }
        return 0;
    }

    void setElement(int i, int j, double elem)
    {
        if (abs(j - i) == 1 || abs(j - i) == 0)
        {
            switch (j - i)
            {
            case -1:
                a[0][j + 1] = elem;
                break;
            case 0:
                a[1][j] = elem;
                break;
            case 1:
                a[2][j - 1] = elem;
                break;
            }
        }
    }

    void print()
    {
        ofstream out("out.txt");
        for (int i = 0; i < a[0].size(); i++)
        {
            for (int j = 0; j < a[0].size(); j++)
                out << getElement(i, j) << " ";
            out << endl;
        }
    }
    int size()
    {
        return a[0].size();
    }
    void clear()
    {
        for (int i = 0; i < a[0].size(); i++)
            for (int j = 0; j < a[0].size(); j++)
                setElement(i, j, 0);
    }
};

class Operations
{
public:
    static void multMatrixByC(vector<vector<double>>& result, const vector<vector<double>>& a, double b)
    {
        for (size_t i{ 0 }; i < a.size(); i++)
            for (size_t j{ 0 }; j < a.size(); j++)
                result[i][j] += a[i][j] * b;
    }

    static void multMatrixByVector(vector<double>& result, Matrix& a, vector<double> b)
    {
        for (size_t i{ 0 }; i < result.size(); i++)
            for (size_t j{ 0 }; j < result.size(); j++)
                result[i] += a.getElement(i, j) * b[j];
    }

    static double norm(const vector<double>& a)
    {
        double result = 0;
        for (int i = 0; i < a.size(); i++)
            result += a[i] * a[i];
        return sqrt(result);
    }
};

class FEM
{
    Matrix a;
    vector<double> q;
    vector<double> qCurrent;
    vector<double> qTime;
    vector<double> b;
    vector<double> grid;
    vector<double> time;
    vector<Element> elements;
    const vector<vector<double>> G
    {
        {1, -1},
        {-1, 1}
    };
    const vector<vector<double>> M
    {
        {2, 1},
        {1, 2}
    };
public:

    double sigma(double x)
    {
        return 1;
    }

    double lambda(double dq, double t, double x)
    {
       return dq;
    }
    
    double difLambda(double dq, double t)
    {
       return 1;
    }

    double lambdaNewton(double dq, double hx, int i, int j, int t)
    {
        double first = difLambda(dq, t);
        double second = (j == 0) ? -1.0 / hx : 1.0 / hx;
        return first * second;
    }

    double dAdQ(double hx, double dq, int i, int j, int r, int t)
    {
        double l1 = lambdaNewton(dq, hx, i, j, t);
        double l2 = lambdaNewton(dq, hx, i, j, t);
        return (i == r) ? (l1 + l2) / (2 * hx) : -(l1 + l2) / (2 * hx);
    }
    
    double f(double x, double t)
    {
       return -8 * x * t * t + x * x;
    }

    double trueQ(double x, double t)
    {
       return x * x * t;
    }
    
    FEM()
    {

    }

    void Init()
    {
        a = Matrix(grid.size());
        b.resize(grid.size());
        q.resize(grid.size());
        qCurrent.resize(grid.size());
        qTime.resize(grid.size());
    }

    void setElements(string fileName)
    {
        ifstream in(fileName);
        int size;
        in >> size;
        //elements.resize(size);
    }

    void setGrid(string fileName, bool isRaw = true)
    {
        ifstream in(fileName);
        double a = 0;
        double b = 1;
        int n = 0;
        in >> a >> b >> n;
        double h = (b - a) / n;
        int i = 0;
        if (isRaw)
        {
           while (true)
           {
              if (a + i * h >= b) break;
              grid.push_back(a + i * h);
              i++;
           }
           grid.push_back(b);
        }
        else
        {
           double pi = 3.1415926535897932;
           for (int i{ 0 }; i < n; i++)
           {
              grid.push_back((b - a) * (1 - cos(pi * i / (2 * n))));
           }
           grid.push_back(b);
        }
        elements.resize(grid.size() - 1);
    }

    void setTime(string fileName, bool isRaw = true)
    {
        ifstream in(fileName);
        double a = 0;
        double b = 0;
        int n = 0;
        in >> a >> b >> n;
        double h = (b - a) / n;
        int i = 0;
        if (isRaw)
        {
           while (true)
           {
              if (a + i * h >= b) break;
              time.push_back(a + i * h);
              i++;
           }
           time.push_back(b);
        }
        else
        {
           double pi = 3.1415926535897932;
           for (int i{ 0 }; i < n; i++)
           {
              time.push_back((b - a)* (1 - cos(pi * i / (2 * n))));
           }
           time.push_back(b);
        }
    }

    void setMatrix(int t)
    {
        for (size_t k{ 0 }; k < elements.size(); k++)
        {
            vector<vector<double>> local(2, vector<double>(2));
            double hk = grid[k + 1] - grid[k];
            double tau = time[t] - time[t - 1];
            double coef{ sigma(grid[k]) * hk / (6.0 * tau) };

            Operations::multMatrixByC(local, M, coef);

            double dq = (qCurrent[k + 1] - qCurrent[k]) / hk;
            coef = (lambda(dq, time[t], grid[k]) + lambda(dq, time[t], grid[k + 1])) / (2 * hk);
            
            Operations::multMatrixByC(local, G, coef);

            a.setElement(k, k, local[0][0] + a.getElement(k, k));
            a.setElement(k, k + 1, local[0][1] + a.getElement(k, k + 1));
            a.setElement(k + 1, k, local[1][0] + a.getElement(k + 1, k));
            a.setElement(k + 1, k + 1, local[1][1] + a.getElement(k + 1, k + 1));


            //set b
            //b[k] += (hk / 6.0) * (f(grid[k],time[t]) * M[0][0] + f(grid[k + 1], time[t]) * M[0][1]);
            //b[k + 1] += (hk / 6.0) * (f(grid[k], time[t]) * M[1][0] + f(grid[k + 1], time[t]) * M[1][1]);

            //test with q

            b[k] += (hk / 6.0) * (f(grid[k], time[t]) * M[0][0] + f(grid[k + 1], time[t]) * M[0][1]) + (hk / (6.0 * tau)) * (qTime[k] * sigma(grid[k]) * M[0][0] + qTime[k + 1] * sigma(grid[k+1]) * M[0][1]);
            b[k + 1] += (hk / 6.0) * (f(grid[k], time[t]) * M[1][0] + f(grid[k + 1], time[t]) * M[1][1]) + (hk / (6.0 * tau)) * (qTime[k] * sigma(grid[k]) * M[1][0] + qTime[k + 1] * sigma(grid[k + 1]) * M[1][1]);
        }
    }

    void setDifMatrix(int t)
    {
        for (size_t k{ 0 }; k < elements.size(); k++)
        {
            vector<vector<double>> local(2, vector<double>(2));
            double hk = grid[k + 1] - grid[k];
            double tau = time[t] - time[t - 1];
            double dq = (qCurrent[k + 1] - qCurrent[k]) / hk;
            double qi = 0;
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    double sum_r = 0;
                    for (int r = 0; r < 2; r++)
                    {
                        if (r == 0)
                        {
                            qi = k;
                        }
                        else
                        {
                            qi = k + 1;
                        }
                        sum_r += dAdQ(hk, dq, i, j, r, t) * qCurrent[qi];
                    }
                    local[i][j] += sum_r;
                    qi = (j == 0) ? k : k + 1;
                    b[k + i] += sum_r * qCurrent[qi];
                }
            }
            a.setElement(k, k, local[0][0] + a.getElement(k, k));
            a.setElement(k, k + 1, local[0][1] + a.getElement(k, k + 1));
            a.setElement(k + 1, k, local[1][0] + a.getElement(k + 1, k));
            a.setElement(k + 1, k + 1, local[1][1] + a.getElement(k + 1, k + 1));
        }
    }

    void setFirstBound(string fileName, int t)
    {
        ifstream in(fileName);
        int elem;
        double value;
        while (in >> elem >> value)
        {
            for (int i = 0; i < a.size(); i++)
            {
                a.setElement(elem - 1, i, 0);
            }
            a.setElement(elem - 1, elem - 1, 1);
            b[elem - 1] = value * time[t];
        }
    }

    void iter(int t, double eps = 1e-6, double w = 1, bool isNewton = false, int maxIterations = 500)
    {
        vector<double> discrepancy(q.size());
        vector<double> discrepancyTrue(q.size());
        vector<double> trQ(q.size());
        double exit = 0;
        double iter = 0;
        do
        {
            a.clear();
            for (int j = 0; j < b.size(); j++)
            {
                q[j] = 0;
                b[j] = 0;
            }
            setMatrix(t);
            if (isNewton) setDifMatrix(t);
            setFirstBound("bound1.txt", t);
            a.print();
            LUMethod lu(a.a, b);
            lu.LU();
            lu.forwardStep();
            lu.backStep();
            lu.getQ(q);
            for (int j = 0; j < q.size(); j++)
                q[j] = w * q[j] + (1 - w) * qCurrent[j];
            for (int j = 0; j < discrepancy.size(); j++)
            {
               trQ[j] = trueQ(grid[j], t);
               discrepancyTrue[j] = abs(q[j] - trQ[j]);
               discrepancy[j] = q[j] - qCurrent[j];
            }
            vector<double> f(q.size());
            Operations::multMatrixByVector(f, a, qCurrent);
            for (int i = 0; i < f.size(); i++)
                f[i] -= b[i];
            exit = Operations::norm(discrepancy) / Operations::norm(q);//Operations::norm(f) / Operations::norm(b);//
            for (int j = 0; j < q.size(); j++)
                qCurrent[j] = q[j];
        } while (abs(exit) > eps && iter++ < maxIterations);
        ofstream out(to_string(t) + ".txt");
        out << iter << endl << 1 - Operations::norm(discrepancyTrue)/Operations::norm(trQ);
        out.close();
    }

    void solve(double eps, bool isNewton, double w, int maxIter)
    {
        for (int i = 0; i < qCurrent.size(); i++)
        {
            qTime[i] = 0;//grid[i] * grid[i];
            qCurrent[i] = qTime[i];
        }
        for (int t = 1; t < time.size(); t++)
        {
            iter(t, eps, w, isNewton, maxIter);
            for (int j = 0; j < q.size(); j++)
                qTime[j] = qCurrent[j];
        }
    }
};

int main()
{
    FEM a;
    bool isRaw = true;
    bool isNewton = false;
    a.setGrid("grid.txt",isRaw);
    a.setTime("time.txt",isRaw);
    a.setElements("elements.txt");
    a.Init();
    a.solve(1e-15, isNewton, 1, 500);
    return 0;
}