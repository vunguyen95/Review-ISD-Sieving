#include"misc.h"
#include"vector.h"
#include<algorithm>
#include<random>
#include<cstdlib>
#include<set>
#include<map>
#include<unordered_map>
#include<chrono>
#include<unordered_set>
using namespace std::chrono;
template < typename SEQUENCE > struct seq_hash //Hash function for vectors
{
	long long int operator() (const SEQUENCE& seq) const
	{
		std::size_t hash = 0;
		//boost::hash_range(hash, seq.begin(), seq.end());
		for (auto &i : seq) { hash ^= i + 0x9e3779b9 + (hash << 6) + (hash >> 2); }
			return hash;
	}
};
using Set = std::unordered_set< vector<int>, seq_hash<vector<int>> >;
using check_Map = std::unordered_map< vector<int>, int, seq_hash<vector<int>> >;
using check_Map1 = std::unordered_map< vector<int>, vector<vector<int>>, seq_hash<vector<int>> >;
typedef map<vector<int>, vector<vector<int>>> Map;
typedef vector<vector<int>> VEC;
/*vector<int> errorGen(int& length, int& weight)
{
	random_device rd;
    	mt19937 gen(rd());
    	uniform_int_distribution<int> distrib(1,length);
    	vector<int> temp;
    	for (int j = 0; j< weight; j++)
    	{
    		auto position = distrib(gen);
    		temp.emplace_back(position);
    	}
	sort(temp.begin(), temp.end());
	//auto error = represent(temp, length);
	return temp;
}*/
vector<int> errorGen(int& length, int& weight)
{
	random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> distrib(1,length);
    vector<int> temp;
    //int position{0};
    for (int j = 0; j< weight; j++)
    {
        auto position = distrib(gen);
        while(find(temp.begin(),temp.end(),position)!= temp.end())
    	{
            position = distrib(gen);
        }
    	temp.emplace_back(position);
	}
	sort(temp.begin(), temp.end());
	return temp;
}
Set samplesGen(Set& res, int& samples, int& length , int& weight)
{
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> distrib(1,length);
    while(res.size()< samples)
    {
        vector<int> temp;
        for(int j = 0; j < weight; j++)
        {
            auto position = distrib(gen);
            temp.emplace_back(position);
        }
        sort(temp.begin(), temp.end());
        auto CheckDuplicate = adjacent_find(temp.begin(), temp.end());
        if(CheckDuplicate == temp.end())
        {
            res.insert(temp);
        }
    }
    return res;
}
//sorting vectors in their first buckets.
void sortSamples(Set& in, map<vector<int>, vector<vector<int>>>& out, int& two_p, int& p)
{
    for(auto sample : in)
    {
        auto getFirst = getComb(sample, two_p, p)[0];
        auto bucketExist = out.find(getFirst);
        if(bucketExist== out.end())
        {
            vector<vector<int>> bucket;
            bucket.emplace_back(sample);
            out.insert(pair< vector<int>, vector<vector<int>> >(getFirst,bucket));
        }
        else
        {
			if(std::find((*bucketExist).second.begin(), (*bucketExist).second.end(), sample) == (*bucketExist).second.end())
            {(*bucketExist).second.emplace_back(sample);
			}
        }
    }
    in.clear();

}
//find next buckets and move vectors. More efficient version

void moveSamples1(Map& m, vector<int>& sample, vector<int>& currentBucket, int p)
{
    vector<int> nextBucket;
    int index1 = p -1;
    int index2 = 1;
    if(currentBucket[0] != sample[p])
    {
        while(index1 >= 0)
        {
            if(currentBucket[index1]!= sample[2*p - index2])
            {
                nextBucket.assign(currentBucket.begin(), currentBucket.begin() + index1);
                auto it = find(sample.begin(), sample.end(), currentBucket[index1]) + 1;
                for(int i = 0; i< index2; i++)
                {
                    nextBucket.emplace_back((*it));
                    it++;
                }
                break;
            }
            else
            {
                index1--;
                index2++;
            }
        }

    }
    if(nextBucket.size()>0)
	{
		auto findWhere = m.find(nextBucket);
		if(findWhere == m.end())
		{
			vector<vector<int>> newBucket;
			newBucket.emplace_back(sample);
			m.insert(pair<vector<int>, vector<vector<int>>>(nextBucket, newBucket));
		}
		else
		{
			if(find((*findWhere).second.begin(), (*findWhere).second.end(), sample) == (*findWhere).second.end()){
			(*findWhere).second.emplace_back(sample);}
		}
	}
}


//combine and move. First look into a bucket, combine if 2 or more found, after combining move the vectors to their corresponding next buckets. in total a vector moves (2p,p) times.
void check(Map& m, VEC& partialParity, vector<int>& partialSyndrome, Set& s, int& round)
{
	vector<int> zero(round);
	for(auto it = m.begin(); it!= m.end(); it++)
	{
		for(int i = 0; i < (*it).second.size(); i++)
		{
			auto syndrome = mul(partialParity, (*it).second[i]);
			if(syndrome == partialSyndrome || syndrome == zero)
			{
				s.insert((*it).second[i]);
			}
		}
	}
}
void combine(Map& m, VEC& partialParity, vector<int>& partialSyndrome, int& two_p, int&p, int& n, Set& s, int& round, int& nbrSamples)
{
	vector<int> zero(round);
	int count{0};
	int duplicate0{0};
    int duplicate1{0};
    Set s0 = s; // Check
    Set s00;
    check_Map1 m00;
    Set s1; // Combine
    Set s11;
    check_Map1 m11;
	while(m.size()>0 && s.size() < nbrSamples)
	{	
		auto it = m.begin();
		auto currentBucket = (*it).first;
		
			for(int i = 0 ; i < (*it).second.size() - 1; i++)
			{
				for(int j = i+1; j < (*it).second.size(); j++)
				{
					auto sum = radd((*it).second[i], (*it).second[j]); //xor samples
					if(sum.size() == two_p)
					{
						auto check = mul(partialParity,sum);
						if(check == partialSyndrome || check == zero)
						{	
							count++;
							s.insert(sum);
							if(s1.find(sum)== s1.end()){
                                        s1.insert(sum);
                                    }
                                    else{
                                        duplicate1++;
                                        s11.insert(sum);
                                        if(m11.find(sum) == m11.end()){
                                            vector<vector<int>> Pairs;
                                            Pairs.emplace_back((*it).second[i]);
                                            Pairs.emplace_back((*it).second[j]);
                                            m11.insert(make_pair(sum, Pairs));
                                        }
                                        else{
                                            auto it1 = m11.find(sum);
                                            (*it1).second.emplace_back((*it).second[i]);
                                            (*it1).second.emplace_back((*it).second[j]);
                                        }
                                    }
                                    //Check duplicates from preliminary check
                                    if(s0.find(sum) != s0.end()){
                                        duplicate0++;
                                        s00.insert(sum);
                                        if(m00.find(sum) == m00.end()){
                                            vector<vector<int>> pairs;
                                            pairs.emplace_back((*it).second[i]);
                                            pairs.emplace_back((*it).second[j]);
                                            m00.insert(make_pair(sum, pairs));
                                        }
                                        else{
                                            auto it2 = m00.find(sum);
                                            (*it2).second.emplace_back((*it).second[i]);
                                            (*it2).second.emplace_back((*it).second[j]);
                                        }
                                    }
						}
					}
				}
			}
		
		for ( auto sample: (*it).second )
		{
			moveSamples1(m, sample, currentBucket, p);
			
		}
		m.erase(it);
		
	}
	cout << "Collision Count: " << log2(count) << endl;
    cout << "Duplicate from preliminary check: " << log2(duplicate0) <<" : " <<duplicate0 << endl;
    cout << "Number of duplicated vectors from check: " << s00.size() << endl;
    cout << "Duplicate from combination: " << log2(duplicate1) <<" : " <<duplicate1 <<endl;
    cout << "Number of duplicated vectors from combination: " << s11.size() << endl;
    vector<vector<int>> max0;
    vector<int> vec0;
    vector<vector<int>> max1;
    vector<int> vec1;
    for(auto it = m00.begin(); it!= m00.end();it++){
        if((*it).second.size() > max0.size()){
            max0 = (*it).second;
            vec0 = (*it).first;
        }
    }
    for(auto it = m11.begin(); it!= m11.end(); it++){
        if((*it).second.size() > max1.size()){
            max1 = (*it).second;
            vec1 = (*it).first;
        }
    }
    cout << "Most duplicated vector from check: " << vec0;
    cout << max0 << endl;
    cout << "Most duplicated vector from combination: " << vec1;
    cout << max1 << endl; 
}



int main()
{
   map<vector<int>, vector<vector<int>>> myMap;
   Set mySet;
   //Parameters:
   int length = 500 + 20 ;// k+ l
   int nbrParities = 20;
   int p = 2;
   int two_p = 4;
   int nbrSamples = pow(2,13.4);
   
   //Parity, random test error e and syndrome:
   auto h0 = parityGen(nbrParities,length);

   vector<vector<int>> H0(h0.size());
		for (int i = 0; i < H0.size(); i++)
		{
			H0[i] = bucketRepresent(h0[i]);
		}


   auto testError = errorGen(length, two_p);
   cout << testError;
   auto syndrome = mul(H0, testError); 
   auto check1 = represent(testError,length);

   //cout << "Parity check H0: " <<endl << h0;
   //cout << "Bucket form: " << endl << H0;
   //cout << "Test error: " << endl << testError;
   cout << "Syndrome:" << endl << syndrome;
   cout << mul1(h0,check1);
   

   
   
   //MERGESET ALGORITHM
   samplesGen(mySet, nbrSamples, length, two_p);
   sortSamples(mySet,myMap,two_p, p);
   
   //auto start = high_resolution_clock::now();

   for(int round = 1; round <= nbrParities; round++)
   {
	cout << "ROUND" <<round <<endl;
   	sortSamples(mySet,myMap, two_p, p);
   	auto partialParity = get(H0, round);
   	vector<int> partialSyndrome;
   	partialSyndrome.assign(syndrome.begin(), syndrome.begin() + round);

	check(myMap, partialParity, partialSyndrome, mySet,round);
	cout << log2(mySet.size()) <<endl;

   	combine(myMap, partialParity, partialSyndrome, two_p, p, length, mySet, round, nbrSamples);
	cout << log2(mySet.size()) << endl;
   }
   
   //auto stop = high_resolution_clock::now();
   //auto duration = duration_cast<microseconds>(stop-start);
   //cout << "TIME :" << duration.count() << endl;

	cout << "DONE" << endl;
	cout << "Check: " << endl;
   auto it = mySet.find(testError);
   if(it!= mySet.end())
   {
	   cout << "FOUND SOLUTION" << endl;
	   cout << testError;
	   cout << *it;
   }
   else
   {cout << "SOLUTION NOT FOUND";}
   
}
