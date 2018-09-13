
def neville(x_data, y_data, x):
    """p = neville(x_data, y_data, x)
    Evaluate the polynomial interpolant p(x) that passes through the
    specified data points by Neville's method."""

    n = len(x_data)
    p = n*[0]
    for k in range(n):
        for i in range(n-k):
            if k == 0:
                p[i] = y_data[i]
            else:
                p[i] = ((x-x_data[i+k])*p[i]+ \
                        (x_data[i]-x)*p[i+1])/ \
                        (x_data[i]-x_data[i+k])
    return p[0]


if __name__ == "__main__":
    x_data = [8.1, 8.3, 8.6, 8.7]
    y_data = [16.9446, 17.56492, 18.50515, 18.82091]
    x = 8.4
    print(neville(x_data, y_data, x))
