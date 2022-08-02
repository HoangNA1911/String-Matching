#include <iostream>
#include <string>
#include <fstream>
#include <cstring>
#include <vector>
using namespace std;

#define d 256
struct wordIF
{
    string name;
    int x;
    int y;
    string dir;
};
void print(char **matrix, int w, int h);
void readFile(char **&matrix, int &width, int &height, vector<string> &word);
vector<wordIF> brute_force(char **matrix, int width, int height, vector<string> &word);
vector<wordIF> rabin_Karp(char **matrix, int width, int height, vector<string> &word, int q);
vector<wordIF> KMP(char **matrix, int width, int height, vector<string> &word);
vector<wordIF> order(char **matrix, int width, int height, vector<string> &word, int q);
void printWPD(vector<wordIF> result);
void LBS(string word, int *&lps);
void readFile(char **&matrix, int &width, int &height, vector<string> &word)
{
    ifstream fin;
    fin.open("input.txt", ios::in);
    string s;
    getline(fin, s, ' ');
    width = stoi(s);
    getline(fin, s);
    height = stoi(s);
    matrix = new char *[height];
    for (int i = 0; i < height; i++)
    {
        matrix[i] = new char[width];
    }
    for (int i = 0; i < height; i++)
    {
        int k = 0;
        // char ch;
        for (int j = 0; j < width - 1; j++)
        {
            getline(fin, s, ' ');
            // strcpy(&ch, s.c_str());
            matrix[i][k] = s[0];
            k++;
        }
        getline(fin, s);
        // strcpy(&ch, s.c_str());
        matrix[i][k] = s[0];
    }
    getline(fin, s);
    while (!fin.eof())
    {
        getline(fin, s);
        if (s == "#")
            break;
        word.push_back(s);
    }
    fin.close();
}
void print(char **matrix, int w, int h)
{
    for (int i = 0; i < h; i++)
    {
        for (int j = 0; j < w; j++)
        {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}
vector<wordIF> brute_force(char **matrix, int width, int height, vector<string> &word)
{
    vector<wordIF> result;
    wordIF dummy;
    int k = 0;
    while (k < word.size())
    {
        string s = word[k];
        int M = s.length();

        for (int i = 0; i < height; i++)
        {
            for (int j = 0; j <= width - M; j++)
            {
                int h;
                for (h = 0; h < M; h++)
                {
                    if (matrix[i][j + h] != s[h])
                        break;
                }

                if (h == M)
                {
                    dummy.x = j;
                    dummy.y = i;
                    dummy.dir = "LR";
                    dummy.name = word[k];
                    result.push_back(dummy);
                }
            }
        }
        for (int j = 0; j < width; j++)
        {
            for (int i = 0; i <= height - M; i++)
            {
                int h;
                for (h = 0; h < M; h++)
                {
                    if (matrix[i + h][j] != s[h])
                        break;
                }
                if (h == M)
                {
                    dummy.x = j;
                    dummy.y = i;
                    dummy.dir = "TD";
                    dummy.name = word[k];
                    result.push_back(dummy);
                }
            }
        }
        k++;
    }
    return result;
}

// ham Rabin karp

vector<wordIF> rabin_Karp(char **matrix, int width, int height, vector<string> &word, int q)
{
    vector<wordIF> result;
    wordIF dummy;
    int k = 0;
    while (k < word.size())
    {
        string s = word[k];
        int M = s.length();
        for (int i = 0; i < height; i++)
        {
            unsigned long long h = 1;
            for (int p = 0; p < M - 1; p++)
                h = (h * d) % q;
            // first hash
            long long a = 0, b = 0;
            for (int p = 0; p < M; p++)
            {
                a = (a * d + s[p]) % q;
                b = (b * d + matrix[i][p]) % q;
            }
            for (int j = 0; j <= width - M; j++)
            {
                int p;
                if (a == b)
                {
                    for (p = 0; p < M; p++)
                    {
                        if (matrix[i][j + p] != s[p])
                            break;
                    }
                    if (p == M)
                    {
                        dummy.x = j;
                        dummy.y = i;
                        dummy.dir = "LR";
                        dummy.name = word[k];
                        result.push_back(dummy);
                    }
                }
                if (j < width - M)
                {
                    b = (d * (b - matrix[i][j] * h) + matrix[i][j + M]) % q;
                    if (b < 0)
                        b = b + q;
                }
            }
        }
        for (int j = 0; j < width; j++)
        {
            unsigned long long h = 1;
            for (int p = 0; p < M - 1; p++)
                h = (h * d) % q;
            // first hash
            long long a = 0, b = 0;
            for (int p = 0; p < M; p++)
            {
                a = (a * d + s[p]) % q;
                b = (b * d + matrix[p][j]) % q;
            }
            for (int i = 0; i <= height - M; i++)
            {
                int p;
                if (a == b)
                {
                    for (p = 0; p < M; p++)
                    {
                        if (matrix[i + p][j] != s[p])
                            break;
                    }
                    if (p == M)
                    {
                        dummy.x = j;
                        dummy.y = i;
                        dummy.dir = "TD";
                        dummy.name = word[k];
                        result.push_back(dummy);
                    }
                }
                if (j < height - M)
                {
                    b = (d * (b - matrix[i][j] * h) + matrix[i][j + M]) % q;
                }
            }
        }
        k++;
    }
    return result;
}

// ham KMP
void LBS(string word, int *&lps)
{
    lps[0] = 0;
    int temp = 0;
    int i = 1;
    while (i < word.size())
    {
        if (word[i] == word[temp])
        {
            temp++;
            lps[i] = temp;
            i++;
        }
        else
        {
            if (temp > 0)
            {
                temp = lps[temp - 1];
            }
            else
            {
                lps[i] = 0;
                i++;
            }
        }
    }
}
vector<wordIF> KMP(char **matrix, int width, int height, vector<string> &word)
{
    vector<wordIF> result;
    wordIF dummy;
    int k = 0;
    while (k < word.size())
    {
        string s = word[k];
        // string s="aaaba";
        int M = s.size();
        int *lps = new int[M];
        LBS(s, lps);
        for (int i = 0; i < height; i++)
        {

            int h = 0; // index cua word
            int j = 0; // index input
            while ((width - j) >= (M - h))
            {
                if (matrix[i][j] == s[h])
                {
                    h++;
                    j++;
                }
                if (h == s.size())
                {
                    dummy.x = j - M;
                    dummy.y = i;
                    dummy.dir = "LR";
                    dummy.name = word[k];
                    result.push_back(dummy);
                    h = lps[h - 1];
                }
                else if (j < width && matrix[i][j] != s[h])
                {
                    if (h != 0)
                        h = lps[h - 1];
                    else
                    {
                        j = j + 1;
                    }
                }
            }
        }

        for (int j = 0; j < width; j++)
        {

            int h = 0; // index cua word
            int i = 0; // index input
            while ((height - i) >= (M - h))
            {
                if (matrix[i][j] == s[h])
                {
                    h++;
                    i++;
                }
                if (h == s.size())
                {
                    dummy.x = j;
                    dummy.y = i - M;
                    dummy.dir = "TD";
                    dummy.name = word[k];
                    result.push_back(dummy);
                    h = lps[h - 1];
                }
                else if (i < height && matrix[i][j] != s[h])
                {
                    if (h != 0)
                        h = lps[h - 1];
                    else
                    {
                        i = i + 1;
                    }
                }
            }
        }
        delete[] lps;
        k++;
    }
    return result;
}
vector<wordIF> order(char **matrix, int width, int height, vector<string> &word, int q)
{
    cout << "1. Brute-force (Naive String-matching) " << endl;
    cout << "2. Rabin-Karp " << endl;
    cout << "3. Knuth-Morris-Pratt " << endl;
    cout << "Nguoi dung nhap thao tac: ";
    int x;
    cin >> x;
    vector<wordIF> checkMatching;
    switch (x)
    {
    case 1:
        checkMatching = brute_force(matrix, width, height, word);
        break;
    case 2:
        checkMatching = rabin_Karp(matrix, width, height, word, q);
        break;
    case 3:
        checkMatching = KMP(matrix, width, height, word);
        break;
    default:
        break;
    }
    return checkMatching;
}
void printWPD(vector<wordIF> result)
{
    cout << result.size() << endl;
    for (int i = 0; i < result.size(); i++)
    {
        cout << result[i].name << " (" << result[i].y << "," << result[i].x << ") " << result[i].dir << endl;
    }
}
int main()
{
    // doc file nguoi dung luu vao mang 2 chieu cap phat dong
    // bao gom: doc so hang va so cot (cap phat dong cho mang 2 chieu va de biet can doc bao nhieu hang trong file de luu vao mang 2 chieu).
    //          doc cac ky tu sau do luu vao mang (luu vao 1 mang 2 chieu kieu char)
    //          doc cac tu can tim kiem trong mang. (luu vao 1 mang 1 chieu string)
    char **matrix;
    int width, height;
    vector<string> word;
    readFile(matrix, width, height, word);
    print(matrix, width, height);
    // tim kiem tren mang, tim kiem theo 3 thuat toan Brute-force, Rabin-Karp, Knuth-Morris-Pratt
    int q = INT_MAX;
    vector<wordIF> checkMatching = order(matrix, width, height, word, q);

    printWPD(checkMatching);

    for (int i = 0; i < height; i++)
    {
        delete[] matrix[i];
    }
    delete[] matrix;
    return 0;
}
