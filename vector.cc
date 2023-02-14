#include"vector.h"
#include<random>
#include<algorithm>

vector<vector<int>> getComb(vector<int>& v, int n, int& r)
{
	vector<int> data(r);
	vector<vector<int>> res;
	combination(v, data, 0, n-1, 0, r, res);
	return res;
}

void combination(vector<int>& v, vector<int>& data, int start, int end, int index, int r, vector<vector<int>>& res)
{
	if(index == r)
	{
		res.push_back(data);
		return;
	}
	for(int i = start; i <= end && end - i + 1 >= r - index; i++)
	{
		data[index] = v[i];
		combination(v, data, i+1, end, index + 1, r, res);
	}
}
vector<char> add( vector<char>& u,  vector<char>& v)//
{
	int n = u.size();
	vector<char> res(n);
	if(u.size() != v.size())
	{
		cout << "Wrong format, check size again";
	}
	else
	{
		for (int i = 0; i < res.size(); i++)
		{
			res[i] = (u[i] ^ v[i]);
		}
	}
	return res;
}
/*vector<int> radd1(vector<int>& v1, vector<int>& v2){
	set<int> s;
	for()
}*/
vector<int> radd(vector<int>& v1, vector<int>& v2)
{
	int i{0},j{0};
	vector<int> res;
	while(i < v1.size() && j < v2.size())
	{
		if(v1[i] >= v2[j] && j == v2.size()-1)
		{
			if (v1[i] > v2[j])
			{
				res.push_back(v2[j]);
			}
			else
			{
				i++;
			}
			while(i < v1.size())
			{
				res.push_back(v1[i]);
				i++;
			}

		}
		else if (v1[i] > v2[j] && j!= v2.size()-1)
		{
			res.push_back(v2[j]);
			j++;
		}
		else if (v1[i] <= v2[j] && i == v1.size()-1)
		{
			if(v1[i] < v2[j])
			{
				res.push_back(v1[i]);
			}
			else
			{
				j++;
			}
			while(j < v2.size())
			{
				res.push_back(v2[j]);
				j++;
			}
		}
		else if (v1[i] < v2[j] && i!= v1.size())
		{
			res.push_back(v1[i]);
			i++;
		}
		else{ i++; j++;}
	}
	return res;
}
int weight(vector<char>& v)
{
	int length = v.size();
	int res = 0;
	for (int i = 0; i < length; i++)
	{
		if (v[i] == 1)
		{
			res += 1;
		}
	}
	return res;
}

vector<int> bucketRepresent(vector<char>& v)
{
	vector<int> res;
	for (int i = 0 ; i< v.size();i++)
	{
		if( v[i]!= 0)
		{res.push_back(i+1);}
	}
	return res;
}
vector<char> represent(vector<int>& v, int& length)
{
	vector<char>res(length);
	for(int i=0; i<v.size();i++)
	{
		res[v[i]-1]=1;
	}
	return res;

}

vector<vector<char>> parityGen(int& l, int&n)
{
	vector<vector<char>> h0(l);
	random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> distrib(0,1);
    for(int i = 0; i< l; i++)
    {
    	h0[i].resize(n);
    	generate(h0[i].begin(), h0[i].end(), [&](){return distrib(gen);});	
    }
    return h0;
    		
	
}
int rscalar(vector<int>& v1, vector<int>& v2)
{
	int res = 0;
	int i{0}, j{0};
	int length = v1.size();
	for (i = 0; i< length; i++)
	{
		while( v2[j] < v1[i])
		{
			if ( j == v2.size()-1)
			{
				return res %2 ;
			}
			j++;
		}
		if(v2[j] == v1[i])
		{
			res++;
		}
	}
	return res % 2;
}
vector<int> mul( vector<vector<int>>& a, vector<int>& b)
{
	int length = a.size();
	vector<int> res(length);
	for (int i = 0 ; i<  res.size(); i++)
	{
		res[i] = rscalar(a[i], b);
	}
	return res;
}

vector<char> mul1( vector<vector<char>>& a, vector<char>& b)// in binary, of the form A.b
{
	int m = a.size();
	int n = a[0].size();
	int p = b.size();
	int i{ 0 }, j{ 0 };
	vector<char> res(m);
	if (n!= p)
		cout << "Wrong format:" << endl;
	else
	{
		while (i < m)
		{
			for (j = 0; j < n; j++)
			{
				res[i] += (a[i][j] * b[j]);
				res[i] = res[i] % 2;
			}
			i++;
		}
	}
	return res;
}
vector<vector<int>> get(vector<vector<int>>& parity, int& i)
{
	vector<vector<int>> res;
	res.assign(parity.begin(), parity.begin() + i);
	return res;
	
}
