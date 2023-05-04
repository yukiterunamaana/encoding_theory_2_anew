import java.util.Arrays;
import java.util.Random;

public class Main {
    public static int[] intToBinary(int n) {
        if (n == 0)
            return new int[] { 0 };
        int[] bits = new int[11];
        int pos = 0;
        while (n > 0) {
            bits[pos++] = n % 2;
            n /= 2;
        }
        int[] result = new int[pos];
        for (int i = 0; i < pos; i++) {
            result[i] = bits[pos - i - 1];
        }
        return result;
    }
    public static int binaryToInt(int[] n) {
        int res = 0;
        for (int i = 0; i < n.length; i++)
            res+=Math.pow(2,n.length-i)*n[i];
        return res;
    }
    public static int[] plus_minus_bin(int[] a, int[] b)
    {
        int[] result = new int[Math.max(a.length, b.length)];
        for (int i = 0; i < result.length; i++) {
            int x = (i < a.length) ? a[i] : 0;
            int y = (i < b.length) ? b[i] : 0;
            result[i] = (x + y) % 2;
        }
        return result;
    }
    public static int[] divide_bin(int[] a, int[] b)
    {
        if (a.length<b.length) //if deg(a)<deg(b)
            return a;
        else
        {
            int diff = a.length-b.length;
            int[] t = new int[diff+1];
            t[0]=1;
            int[] temp = b;
            System.out.println(Arrays.toString(t));
            temp=multiply_bin(t,temp);
            System.out.println(Arrays.toString(temp));
            return divide_bin(temp, b);
        }
    }
    public static int[] multiply_bin(int[] a, int[] b)
    {
        int[] result = new int[a.length + b.length - 1];
        for (int i = 0; i < a.length; i++)
            for (int j = 0; j < b.length; j++)
                result[i + j] ++;
        //must be [1,2,1,0], why the extra power?!
        for (int k=0; k<result.length; k++)
            result[k]%=2;
//        int i = 0;
//        while (i < result.length && result[i] == 0) i++;
//        int[] trimmed = new int[i+1];
//        System.arraycopy(result, 0, trimmed, 0, i);
        return result;//divide_bin(trimmed,intToBinary(Q_primal));
    }

    static final int Q = 1024;
    static final int Q_primal = 1033;
    static final int k=5;
    static final int n=7;
    static final int d = n-k+1;
    static final int t = (int) Math.floor(((double)d-1)/2);
    static int calcPolynomeValueANEW(int x, int[] p)
    {
        int res=0;
        for (int i=0; i<p.length; i++)
            res+=p[p.length-1-i]*Math.pow(x,p.length-1-i);
        return res % Q;
    }
    static int[] applyErrors(int[] p)
    {
        Random r = new Random();
        for (int i = 0; i<t; i++)
        {
            int j = Math.abs(r.nextInt() % p.length-1);
            p[j] = Math.abs(r.nextInt() % Q);
        }
        return p;
    }

    static int[] gaoEncode(int[] m, int[] a)
    {
        int[] c = new int[a.length];
        for (int i=0;i<a.length;i++)
            c[i]=calcPolynomeValueANEW(a[i],m);
        applyErrors(c);
        return c;
    }

    //    public static int[] polynominalG0(int[] a) {
//        int n = a.length;
//        int[] poly = new int[n+1];
//        poly[0] = 1;
//        for (int i = 0; i < n; i++) {
//            int[] prevPoly = poly.clone();
//            poly[0] = gf.subtract(0, a[i]); // x - a[i]
//            for (int j = 1; j <= i + 1; j++) {
//                poly[j] = gf.add(gf.multiply(prevPoly[j-1], poly[0]), prevPoly[j]);
//            }
//        }
//        return poly;
//    }
    public static void main(String[] args) {
        System.out.printf("d = %d, t = %d\n",d,t);

        int[] m = {2,10,3,5,1}; //x4 + 5*x3 + 3*x2 + 10*x + 2
        int[] a = {1,2,3,4,5,6,7};
        //System.out.println(Arrays.toString(gaoEncode(m,a)));

        //System.out.println(Arrays.toString(intToBinary(Q_primal)));
        int[] binary_ten=multiply_bin(intToBinary(6),intToBinary(3)); //should be [1,0,1,0]
        System.out.println(Arrays.toString(binary_ten));
        System.out.println(binaryToInt(binary_ten));




    }
}
