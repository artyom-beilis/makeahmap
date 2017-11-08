__kernel void test(int N,__global float4 const *A,__global float4 const *B,__global float4 *C)
{
    int r = get_global_id(0);
    int c = get_global_id(1);
    //if(r < N && c < N) {
        int pos = r*(N/4)+c;
        C[pos] = tanh(A[pos]*A[pos] + B[pos]*B[pos]);
    //}
}
