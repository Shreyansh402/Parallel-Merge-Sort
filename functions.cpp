// Assignment - 1 COL730
// Shreyansh Jain
// 2021MT10230

#include "functions.h"
#include <iostream>
#include <omp.h>
#include <vector>

using namespace std;

void parallel_merge_sort_helper(vector<Record *> &arr, int start, int end);
void merge(vector<Record *> &arr, int start, int mid, int end);
void merge_sort_helper(vector<Record *> &arr, int start, int end);
int parallel_binary_search_helper(vector<Record> &arr, int key, int start, int end);
int binary_search_helper(vector<Record> &arr, int key, int start, int end);

const int THRESHOLD = 1024;

void merge_sort(vector<Record> &arr)
{
    if (arr.size() < 2)
        return;
    vector<Record *> arr_ptr(arr.size());
    for (int i = 0; i < arr.size(); i++)
    {
        arr_ptr[i] = &arr[i];
    }
    merge_sort_helper(arr_ptr, 0, arr.size());
    vector<Record> temp(arr.size());
    for (int i = 0; i < arr.size(); i++)
    {
        temp[i] = *arr_ptr[i];
    }
    arr = temp;
}
void merge_sort_helper(vector<Record *> &arr, int start, int end)
{
    if (end - start < 2)
        return;
    int mid = (start + end) / 2;
    merge_sort_helper(arr, start, mid);
    merge_sort_helper(arr, mid, end);
    merge(arr, start, mid, end);
}

void parallel_merge_sort(vector<Record> &arr)
{
    // create a vector of pointers to the records
    vector<Record *> arr_ptr(arr.size());
    vector<Record> temp(arr.size());
    for (int i = 0; i < arr.size(); i++)
    {
        arr_ptr[i] = &arr[i];
    }
#pragma omp parallel default(none) shared(arr_ptr, temp, arr)
    {
#pragma omp single
        parallel_merge_sort_helper(arr_ptr, 0, arr_ptr.size());
        // copy the sorted records back to the original vector

#pragma omp for schedule(static)
        for (int i = 0; i < arr.size(); i++)
        {
            temp[i] = *arr_ptr[i];
        }
        arr = temp;
    }
}

void parallel_merge_sort_helper(vector<Record *> &arr, int start, int end)
{
    if (end - start < 2)
        return;
    if (end - start < THRESHOLD) // Sequential sort for small arrays
    {
        merge_sort_helper(arr, start, end);
        return;
    }
    int mid = (start + end) / 2;
// #pragma omp parallel default(none) shared(arr, mid, left, right)
#pragma omp task shared(arr) // Add shared(arr) is important
    parallel_merge_sort_helper(arr, start, mid);
#pragma omp task shared(arr)
    parallel_merge_sort_helper(arr, mid, end);
#pragma omp taskwait
    merge(arr, start, mid, end);
}

void merge(vector<Record *> &arr, int start, int mid, int end)
{
    // vector<Record> temp(arr.begin() + start, arr.begin() + end);
    vector<Record *> temp(end - start);
    int i = start, j = mid, k = 0;
    while (i < mid && j < end)
    {
        if (arr[i]->key < arr[j]->key)
            temp[k++] = arr[i++];
        else
            temp[k++] = arr[j++];
    }
    while (i < mid)
        temp[k++] = arr[i++];
    while (j < end)
        temp[k++] = arr[j++];
    // for (int i = start; i < end; i++)
    // {
    //     arr[i] = temp[i - start];
    // }
    std::copy(temp.begin(), temp.end(), arr.begin() + start);
}

vector<Record> binary_search(vector<Record> &arr, int key)
{
    vector<Record> result;
    int left = binary_search_helper(arr, key, 0, arr.size() - 1);
    // printf("Left: %d\n", left);
    if (left != -1)
    {
        int right = binary_search_helper(arr, key + 1, 0, arr.size() - 1); // using 0 as start
        if (right == -1)
            right = arr.size();
        // printf("Right: %d\n", right);
        result.assign(arr.begin() + left, arr.begin() + right);
    }
    return result;
}

int binary_search_helper(vector<Record> &arr, int key, int start, int end)
{
    if (start > end || arr[start].key > key || arr[end].key < key)
        return -1;
    if (arr[start].key == key)
        return start;
    int mid;
    while (start <= end)
    {
        mid = (start + end) / 2;
        if (arr[mid].key < key)
        {
            start = mid + 1;
        }
        else
        {
            end = mid - 1;
        }
    }
    return start;
}

// vector<Record> parallel_binary_search(vector<Record> &arr, int key)
// {
//     vector<Record> result;
//     int left, right;
// #pragma omp parallel default(none) shared(arr, result, left, right) firstprivate(key)
//     {
// #pragma omp task
//         {
//             left = parallel_binary_search_helper(arr, key, 0, arr.size() - 1);
//             // printf("Left: %d\n", left);
//         }
// #pragma omp task
//         {
//             right = parallel_binary_search_helper(arr, key + 1, 0, arr.size() - 1);
//             // printf("Right: %d\n", right);
//             if (right == -1)
//                 right = arr.size();
//         }
//     }
//     if (left == -1)
//     {
//         return result;
//     }
//     result.assign(arr.begin() + left, arr.begin() + right);
//     return result;
// }

vector<Record> parallel_binary_search(vector<Record> &arr, int key)
{
    vector<Record> result;
    int left = arr.size();
    int right = arr.size();
#pragma omp parallel default(none) shared(arr, result, left, right) firstprivate(key)
    {
        int nthreads = omp_get_num_threads();
        int i = omp_get_thread_num();
        // divide the array into nthreads parts and search in parallel
        int new_start = i * arr.size() / nthreads;
        int new_end = (i + 1) * arr.size() / nthreads - 1;
        int lefter = binary_search_helper(arr, key, new_start, new_end);
        if (lefter != -1)
        {
            lefter = min(left, lefter);
#pragma omp atomic write
            left = lefter;
        }
        // divide the array into nthreads parts and search in parallel
        int righter = binary_search_helper(arr, key + 1, new_start, new_end);
        if (righter != -1)
        {
            righter = min(right, righter);
#pragma omp atomic write
            right = righter;
        }
    }

    if (left == arr.size())
    {
        return result;
    }

    result.assign(arr.begin() + left, arr.begin() + right);
    return result;
}

// int parallel_binary_search_helper(vector<Record> &arr, int key, int start, int end)
// {
//     if (start > end || arr[start].key > key || arr[end].key < key)
//         return -1;
//     if (arr[start].key == key)
//         return start;
//     if (end - start < 100000)
//     {
//         return binary_search_helper(arr, key, start, end);
//     }
//     int mid = (start + end) / 2;
//     int nthreads = omp_get_num_threads();
//     // int nthreads = 8;
//     int ans = end + 1;
//     if (arr[mid].key < key)
//     {
//         for (int i = 0; i < nthreads; i++)
//         {
// #pragma omp task default(none) shared(ans, arr, key) firstprivate(i, nthreads, end, mid)
//             {
//                 int new_start = i * (end - mid) / nthreads + mid;
//                 int new_end = (i + 1) * (end - mid) / nthreads + mid - 1;
//                 int left = parallel_binary_search_helper(arr, key, new_start, new_end);
//                 if (left != -1)
//                 {
//                     left = min(left, ans);
// #pragma omp atomic write
//                     ans = left;
//                 }
//             }
//         }
//     }
//     else
//     {
//         for (int i = 0; i < nthreads; i++)
//         {
// #pragma omp task default(none) shared(ans, arr, key) firstprivate(i, nthreads, start, mid)
//             {
//                 int new_start = i * (mid - start) / nthreads + start;
//                 int new_end = (i + 1) * (mid - start) / nthreads + start - 1;
//                 int left = parallel_binary_search_helper(arr, key, new_start, new_end);
//                 if (left != -1)
//                 {
//                     left = min(left, ans);
// #pragma omp atomic write
//                     ans = left;
//                 }
//             }
//         }
//     }
// #pragma omp taskwait
//     if (ans == end + 1)
//         return -1;
//     return ans;
// }

int main()
{
    return 0;
}