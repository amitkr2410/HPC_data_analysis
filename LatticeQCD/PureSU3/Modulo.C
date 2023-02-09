int Modulo(int a, int b) 
{
    int m = a % b;
    if (m < 0) {
      // m += (b < 0) ? -b : b; // avoid this form: it is UB when b == INT_MIN
      m = (b < 0) ? m - b : m + b;
    }
    return m;
}
