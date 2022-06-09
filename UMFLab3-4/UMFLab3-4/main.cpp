#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <set>

using namespace std;

//Lab4

//class Matrix
//{
//   //vector<vector<double>> a;
//   vector<double> a;
//   int n;
//public:
//
//   Matrix()
//   {
//
//   }
//
//   Matrix(Matrix& mat)
//   {
//      a.resize(mat.n * mat.n);
//      n = mat.n;
//      for (int i = 0; i < a.size(); i++)
//      {
//         a[i] = mat.a[i];
//      }
//      /* a.resize(mat.n * mat.n);
//       for (int i = 0; i < a.size(); i++)
//       {
//          a[i].resize(mat.size(i));
//          for (int j = 0; j < a[i].size(); j++)
//          {
//             a[i][j] = mat.getElement(i, j);
//          }
//       }*/
//   }
//
//   Matrix(int n)
//   {
//      a.resize(n * n);
//      this->n = n;
//      //a.resize(n, vector<double>(n));
//   }
//
//   double getElement(int i, int j)
//   {
//      return a[j * n + i];//a[i][j];
//   }
//
//   void setElement(int i, int j, double value)
//   {
//      //a[i][j] = value;
//      a[j * n + i] = value;
//   }
//
//   void print()
//   {
//      ofstream out("out.txt");
//      for (int i = 0; i < n; i++)
//      {
//         for (int j = 0; j < n; j++)
//         {
//            out << setprecision(5) << fixed << getElement(i, j) << " ";
//         }
//         out << endl;
//      }
//   }
//
//   int size()
//   {
//      return n;
//   }
//   int size(int i)
//   {
//      return n;
//   }
//};

class Matrix
{
   vector<int> ig;
   vector<int> jg;
   vector<double> gu;
   vector<double> gl;
   vector<double> di;
public:
   Matrix()
   {

   }
   Matrix(Matrix& mat)
   {
      jg = vector<int>(mat.jg);
      ig = vector<int>(mat.ig);
      di = vector<double>(mat.di);
      gu = vector<double>(mat.gu);
      gl = vector<double>(mat.gl);
   }

   Matrix(int n)
   {

   }

   Matrix(vector<int>& i, vector<int>& j, int n)
   {
      ig = vector<int>(i);
      jg = vector<int>(j);
      gu.resize(ig[2 * n]);
      gl.resize(ig[2 * n]);
      di.resize(2 * n);
   }

   double getElement(int i, int j)
   {
      if (i == j) return di[i];
      
      if (i > j)
      {
         for (int ind = ig[i]; ind < ig[i + 1]; ind++)
         {
            if (jg[ind] == j)
            {
               return gl[ind];
            }
         }
      }
      
      if (i < j)
      {
         for (int ind = ig[j]; ind < ig[j + 1]; ind++)
         {
            if (jg[ind] == i)
            {
               return gu[ind];
            }
         }
      }
      return 0;
   }
   void clear()
   {
      int size = di.size();
      di.clear();
      di.resize(size);
      size = gu.size();
      gu.clear();
      gu.resize(size);
      size = gl.size();
      gl.clear();
      gl.resize(size);
   }
   void Tranc()
   {
      auto tmp = vector<double>(gl);
      gl = vector<double>(gu);
      gu = vector <double>(tmp);
   }
   void setElement(int i, int j, double value)
   {
      if (i == j) di[i] = value;
      if (i > j)
      {
         for (int ind = ig[i]; ind < ig[i + 1]; ind++)
         {
            if (jg[ind] == j)
            {
               gl[ind] = value;
            }
         }
      }
      if (i < j)
      {
         for (int ind = ig[j]; ind < ig[j + 1]; ind++)
         {
            if (jg[ind] == i)
            {
               gu[ind] = value;
            }
         }
      }
   }

   vector<int> getColnsArray()
   {
      return jg;
   }

   vector<int> getRowsArray()
   {
      return ig;
   }

   int size()
   {
      return di.size();
   }

   int size(int i)
   {
      return di.size();
   }

   void print(string filename = "out.txt")
   {
      ofstream out(filename);
      for (int i = 0; i < di.size(); i++)
      {
         for (int j = 0; j < di.size(); j++)
         {
            out << setprecision(5) << fixed << getElement(i, j) << " ";
         }
         out << endl;
      }
   }
};

class BCG
{
   Matrix a;
   //vector<vector<double>> a;
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

   void multMatrixByVector(vector<double>& res, Matrix& mat, vector<double>& v)
   {
      auto jg = mat.getColnsArray();
      auto ig = mat.getRowsArray();

      for (int i = 0; i < res.size(); i++)
      {
         for (int indI = ig[i]; indI < ig[i + 1] - 2; indI++)
         {
            int indJ = jg[indI];
            res[i] += mat.getElement(i, indJ) * v[indJ];
            res[indJ] += mat.getElement(indJ, i) * v[i];
         }
         res[i] += mat.getElement(i, i) * v[i];
      }
   }

   double scalar(vector<double>& a, vector<double>& b)
   {
      double sum = 0;
      for (int i = 0; i < a.size(); i++)
         sum += a[i] * b[i];
      return sum;
   }

   Matrix aT()
   {
      Matrix res(a);
      res.Tranc();
      return res;
   }
public:

   BCG(Matrix& mat, vector<double>& f, int iterations)
   {
      a = Matrix(mat);
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
   double sigma;
   double hi;
   double lambda;
   double omega;

   Element()
   {
      number = 0;
      nodes.resize(4);
      sigma = hi = lambda = omega = 0;
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

class LOS
{
   Matrix a;
   vector<double> x;
   vector<double> b;
   vector<double> r;
   vector<double> z;
   vector<double> p;
   int maxIterations;

   double norm(vector<double>& x)
   {
      double res = 0;
      for (int i = 0; i < x.size(); i++)
         res += x[i] * x[i];
      return sqrt(res);
   }

   void multMatrixByVector(vector<double>& res, Matrix& mat, vector<double>& v)
   {
      auto jg = mat.getColnsArray();
      auto ig = mat.getRowsArray();

      for (int i = 0; i < res.size(); i++)
      {
         for (int indI = ig[i]; indI < ig[i + 1] - 2; indI++)
         {
            int indJ = jg[indI];
            res[i] += mat.getElement(i, indJ) * v[indJ];
            res[indJ] += mat.getElement(indJ, i) * v[i];
         }
         res[i] += mat.getElement(i, i) * v[i];
      }
   }

   double scalar(vector<double>& a, vector<double>& b)
   {
      double sum = 0;
      for (int i = 0; i < a.size(); i++)
         sum += a[i] * b[i];
      return sum;
   }
public:
   LOS()
   {

   }
   LOS(Matrix& mat, vector<double>& b, int iter)
   {
      a = Matrix(mat);
      this->b = vector<double>(b);
      r.resize(b.size());
      z.resize(b.size());
      p.resize(b.size());
      x.resize(b.size());
      maxIterations = iter;
   }

   void calculate(vector<double>& x_start, double eps = 1e-8)
   {
      multMatrixByVector(r, a, x_start);
      for (int i = 0; i < r.size(); i++)
      {
         r[i] = b[i] - r[i];
         z[i] = r[i];
      }
      multMatrixByVector(p, a, z);
      double disp = 1;
      for (int k = 1; k <= maxIterations && disp > eps; k++)
      {
         double alpha = scalar(p, r) / scalar(p, p);
         for (int i = 0; i < x.size(); i++)
         {
            x[i] += alpha * z[i];
            r[i] -= alpha * p[i];
         }
         vector<double> tmp(b.size());
         multMatrixByVector(tmp, a, r);
         double beta = scalar(p, tmp) / scalar(p, p);
         for (int i = 0; i < z.size(); i++)
         {
            z[i] = r[i] + beta * z[i];
            p[i] = tmp[i] + beta * p[i];
         }
         disp = norm(r);
      }
   }
   void getResult(vector<double>& res)
   {
      for (int i = 0; i < x.size(); i++)
         res[i] = x[i];
   }
};

class LU
{
   Matrix lu;
   vector<double> b;
   vector<double> x;
   vector<double> y;

   void forward()
   {
      auto jg = lu.getColnsArray();
      auto ig = lu.getRowsArray();

      for (int i = 0; i < lu.size(); i++)
      {
         double sum = 0;
         for (int indI = ig[i]; indI < ig[i + 1] - 2; indI++)
         {
            int indJ = jg[indI];
            if (indJ >= i) break;
            sum += y[indJ] * lu.getElement(i, indJ);
         }
         y[i] = b[i] - sum;
      }
   }

   void backward()
   {
      auto jg = lu.getColnsArray();
      auto ig = lu.getRowsArray();

      for (int i = lu.size() - 1; i >= 0; i--)
      {
         double sum = 0;
         //for (int indI = ig[i]; indI < ig[i + 1] - 2; indI++)
         for (int indI = ig[i + 1] - 3; indI >= ig[i]; indI--)
         {
            int indJ = jg[indI];
            if (indJ <= i) break;
            sum += x[indJ] * lu.getElement(i, indJ);
         }
         x[i] = (y[i] - sum) / lu.getElement(i, i);
      }
   }
public:
   LU(Matrix& a, vector<double>& b)
   {
      lu = Matrix(a);
      this->b = vector<double>(b);
      x.resize(b.size());
      y.resize(b.size());
      auto jg = a.getColnsArray();
      auto ig = a.getRowsArray();

      for (int i = 0; i < a.size(); i++)
      {
         for (int indI = ig[i]; indI < ig[i + 1] - 2; indI++)
         {
            int indJ = jg[indI];
            if (i <= indJ)
            {
               double sum = 0;
               for (int k = 0; k < i; k++)
                  sum += lu.getElement(i, k) * lu.getElement(k, indJ);
               lu.setElement(i, indJ, a.getElement(i, indJ) - sum);
            }
            if (i > indJ)
            {
               double sum = 0;
               for (int k = 0; k < indJ; k++)
               {
                  double ukj = lu.getElement(k, indJ);
                  sum += lu.getElement(i, k) * ukj;
               }
               lu.setElement(i, indJ, (a.getElement(i, indJ) - sum) / lu.getElement(indJ, indJ));
            }
         }
         
      }
   }
   
   void print(string filename)
   {
      lu.print(filename);
   }

   void solve()
   {
      forward();
      backward();
   }
   void getResult(vector<double>& res)
   {
      for (int i = 0; i < x.size(); i++)
         res[i] = x[i];
   }
};



class FEM
{
   //vector<vector<Point>> grid;
   vector<double> x;
   vector<double> y;
   vector<Element> elements;
   Matrix a;
   vector<double> b;
   vector<double> q;
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

   int getGlobalNum(int i, int j, int k, int nx) // i - global num from mesh, j - local num, nx - count of elements in row
   {
      return (i + k) * nx + j;
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
      a.setElement(2 * gi, 2 * gj, a.getElement(2 * gi, 2 * gj) + local[0][0]);
      a.setElement(2 * gi, 2 * gj + 1, a.getElement(2 * gi, 2 * gj + 1) + local[0][1]);
      a.setElement(2 * gi + 1, 2 * gj, a.getElement(2 * gi + 1, 2 * gj) + local[1][0]);
      a.setElement(2 * gi + 1, 2 * gj + 1, a.getElement(2 * gi + 1, 2 * gj + 1) + local[1][1]);
      b[2 * gi] += blocal[0];
      b[2 * gi + 1] += blocal[1];
   }

   Point getPoint(int i, int j)
   {
      return Point(x[j], y[i]);
   }

   void makePortrait()
   {
      vector<set<int>> list(x.size() * y.size());

      for (int k = 0; k < elements.size(); k++)
      {
         for (int i = 0; i < 4; i++)
         {
            for (int j = 0; j < 4; j++)
            {
               int ind1 = elements[k].getNode(i);
               int ind2 = elements[k].getNode(j);
               if (ind1 < ind2) swap(ind1, ind2);
               list[ind1].insert(ind2);
            }
         }
      }
      vector<int> ig(2 * x.size() * y.size() + 1);
      ig[0] = ig[1] = 0;
      ig[2] = 1;
      for (int i = 1; i < x.size() * y.size(); i++)
      {
         ig[2 * i + 1] = ig[2 * i] + list[i].size() * 2;
         ig[2 * (i + 1)] = ig[2 * i + 1] + list[i].size() * 2 + 1;
      }
      vector<int> jg(ig[2 * x.size() * y.size()]);
      for (int i = 1, k = 1; i < x.size() * y.size(); i++)
      {
         for (int j : list[i])
         {

            jg[k] = 2 * j;
            jg[k + 1] = 2 * j + 1;
            k += 2;
         }
         for (int j : list[i])
         {

            jg[k] = 2 * j;
            jg[k + 1] = 2 * j + 1;
            k += 2;
         }
         jg[k] = 2 * i;
         k++;
      }
      a = Matrix(ig, jg, x.size() * y.size());
      //a = Matrix(2 * y.size() * x.size());
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
      for (int i = 0; i < ny; i++)
         y.push_back(ay + hy * i);
      y.push_back(by);
      for (int j = 0; j < nx; j++)
         x.push_back(ax + hx * j);
      x.push_back(bx);

      //This method need a lots of memory, Now using 2 vectors for x and y points + func getPoint by i and j
      /* for (int i = 0; i <= ny; i++)
       {
          vector<Point> tmp;
          for (int j = 0; j < nx; j++)
          {
             tmp.push_back(Point(ax + hx * j, ay + hy * i));

          }
          tmp.push_back(Point(bx, ay + hy * i));
          grid.push_back(tmp);
       }*/
      int num = 0;
      for (int i = 0; i < y.size() - 1; i++)
         for (int j = 0; j < x.size() - 1; j++)
         {
            elements[num].number = num;
            elements[num++].setNodes(i, nx);
         }
      b.resize(2 * y.size() * x.size());
      q.resize(2 * y.size() * x.size());
      makePortrait();
   }

   void setParams(string filename = "", double sig = 0, double omeg = 0, double khi = 0, double lam = 0)
   {
      if (filename == "")
      {
         for (auto& a : elements)
         {
            a.sigma = sig;
            a.omega = omeg;
            a.hi = khi;
            a.lambda = lam;
         }
      }
      else
      {
         ifstream in(filename);
         if (in.is_open())
         {
            int k;
            in >> k;
            if (k == elements.size())
            {
               double sigma, omega, hi, lambda;
               int i = 0;
               while (in >> sigma >> omega >> hi >> lambda)
               {
                  elements[i].sigma = sigma;
                  elements[i].omega = omega;
                  elements[i].hi = hi;
                  elements[i].lambda = lambda;
               }
            }
         }
      }
   }

   void setMatrix()
   {
      vector<vector<double>> aij(2, vector<double>(2));
      vector<double> bij(2);
      vector<vector<double>> pij(M.size(), vector<double>(M.size()));
      vector<vector<double>> tmp(M.size(), vector<double>(M.size()));
      vector<vector<double>> cij(M.size(), vector<double>(M.size()));
      vector<vector<double>> vec(M.size(), vector<double>(M.size()));
      int nx = x.size();//grid[0].size();
      for (int k = 0; k < elements.size(); k++)
      {
         double hx = abs(getPoint(elements[k].getNode(0) / nx, getLocalNum(0, nx, elements[k].getNode(0))).X - getPoint(elements[k].getNode(1) / nx, getLocalNum(1, nx, elements[k].getNode(1))).X);//abs(grid[elements[k].getNode(0) / nx][getLocalNum(0, nx, elements[k].getNode(0))].X - grid[elements[k].getNode(1) / nx][getLocalNum(1, nx, elements[k].getNode(1))].X);
         double hy = abs(getPoint(elements[k].getNode(0) / nx, getLocalNum(0, nx, elements[k].getNode(0))).Y - getPoint(elements[k].getNode(2) / nx, getLocalNum(2, nx, elements[k].getNode(2))).Y);//abs(grid[elements[k].getNode(0) / nx][getLocalNum(0, nx, elements[k].getNode(0))].Y - grid[elements[k].getNode(2) / nx][getLocalNum(2, nx, elements[k].getNode(2))].Y);
         double coef = elements[k].omega * elements[k].sigma * hx * hy / 36.0;

         multMatrixByConst(cij, M, coef);
         multMatrixByConst(vec, M, hx * hy / 36.0);
         multMatrixByConst(pij, G, elements[k].lambda * hx / (6 * hy));
         multMatrixByConst(tmp, G, elements[k].lambda * hy / (6 * hx));
         AddMatrixs(pij, pij, tmp);
         multMatrixByConst(tmp, M, elements[k].hi * elements[k].omega * elements[k].omega);
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
               double x = getPoint(li, lj).X;
               double y = getPoint(li, lj).Y;
               bij[0] = fs(x, y) * vec[i][j];
               bij[1] = fc(x, y) * vec[i][j];
               addLocalInGlobal(gi, gj, aij, bij);
            }
         }
      }
      //a.print();
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
      return -1;
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
            for (int i = 0; i < x.size(); i++)
            {
               for (int j = 0; j < a.size(i); j++)
               {
                  a.setElement(2 * i, j, 0);
                  a.setElement(2 * i + 1, j, 0);
               }
               a.setElement(2 * i, 2 * i, 1);
               a.setElement(2 * i + 1, 2 * i + 1, 1);
               b[2 * i] = GetBound(getPoint(0, i).X, getPoint(0, i).Y, num, 0);
               b[2 * i + 1] = GetBound(getPoint(0, i).X, getPoint(0, i).Y, num, 1);

            }
            break;
         case 2:
            for (int i = 0; i < y.size() * x.size(); i += x.size())
            {
               for (int j = 0; j < a.size(i); j++)
               {
                  a.setElement(2 * i, j, 0);
                  a.setElement(2 * i + 1, j, 0);
               }
               a.setElement(2 * i, 2 * i, 1);
               a.setElement(2 * i + 1, 2 * i + 1, 1);
               int q = i / x.size();
               b[2 * i] = GetBound(getPoint(q, 0).X, getPoint(q, 0).Y, num, 0);
               b[2 * i + 1] = GetBound(getPoint(q, 0).X, getPoint(q, 0).Y, num, 1);
            }
            break;
         case 3:
            for (int k = 0; k < x.size(); k++)
            {
               int i = (y.size() - 1) * x.size() + k;
               for (int j = 0; j < a.size(i); j++)
               {
                  a.setElement(2 * i, j, 0);
                  a.setElement(2 * i + 1, j, 0);
               }
               a.setElement(2 * i, 2 * i, 1);
               a.setElement(2 * i + 1, 2 * i + 1, 1);
               b[2 * i] = GetBound(getPoint(y.size() - 1, k).X, getPoint(y.size() - 1, k).Y, num, 0);
               b[2 * i + 1] = GetBound(getPoint(y.size() - 1, k).X, getPoint(y.size() - 1, k).Y, num, 1);
            }
            break;
         case 4:
            for (int i = x.size() - 1; i < y.size() * x.size() + x.size() - 1; i += x.size())
            {
               for (int j = 0; j < a.size(i); j++)
               {
                  a.setElement(2 * i, j, 0);
                  a.setElement(2 * i + 1, j, 0);
               }
               a.setElement(2 * i, 2 * i, 1);
               a.setElement(2 * i + 1, 2 * i + 1, 1);
               int q = i / x.size();
               b[2 * i] = GetBound(getPoint(q, x.size() - 1).X, getPoint(q, x.size() - 1).Y, num, 0);
               b[2 * i + 1] = GetBound(getPoint(q, x.size() - 1).X, getPoint(q, x.size() - 1).Y, num, 1);
            }
            break;
         }
      }
      //a.print();
   }

   void solve()
   {
      //BCG solv(a, b, 2);
      LU testi(a, b);
      a.print();
      testi.print("lu.txt");
      testi.solve();
      testi.getResult(q);
      /*LOS solv(a, b, 10000);
      auto x = vector<double>(b.size(), 0);
      solv.calculate(x);
      solv.getResult(q);*/
   }

   void print()
   {
      ofstream out("solution.txt");
      for (int i = 0; i < y.size(); i++)
         for (int j = 0; j < x.size(); j++)
         {
            int k = i * x.size() + j;
            out << q[2 * k] << endl;
         }

   }
};

int main()
{
   FEM fem;
   fem.Init(2, 2);
   fem.setParams("", 0, 0, 0, 1);
   fem.setMatrix();
   fem.FirstBound();
   fem.solve();
   fem.print();
}