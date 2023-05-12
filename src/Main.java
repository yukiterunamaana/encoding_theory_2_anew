import java.util.*;

class Commons
{
    static int[] trim(int[] arr) //V
    {
        int c = 0;
        for (int i = arr.length-1; i>=0; i--)
            if (arr[i]==0)
                c++;
            else
                break;
        int[] res = Arrays.copyOfRange(arr,0,arr.length - c);
        if (res.length>0)
            return res;
        else return new int[] {0};
    }

}

class Task1
{
    static final int Q = 8;
    static final int Q_primal = 11;


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
            result[i] = bits[i];
        }
        return result;
    } //V
    public static int binaryToInt(int[] n) {
        int res = 0;
        for (int i = 0; i < n.length; i++)
            res+=Math.pow(2,i)*n[i];
        return res;
    } //V

    static int[] plus_minus_bin_(int[] a, int[] b) //V
    {
        int[] result = new int[Math.max(a.length, b.length)];
        for (int i = 0; i < result.length; i++) {
            int x = (i < a.length) ? a[i] : 0;
            int y = (i < b.length) ? b[i] : 0;
            result[i] = (x + y) % 2;
        }
        return Commons.trim(result);
    }
    static int[][] divide_bin_(int[] a, int[] b, int[] q) //V
    {
        int aDeg = a.length - 1;
        int bDeg = b.length - 1;
        if (aDeg < bDeg) return new int[][]{q,a};
        int[] quotient = new int[aDeg - bDeg + 1];
        quotient[quotient.length-1]=1;
        q=plus_minus_bin_(q,quotient);
        int[] sub = temp_mult_bin(b,quotient);
        int[] temp = plus_minus_bin_(a, sub);
        temp=Commons.trim(temp);
        return divide_bin_(temp, b, q);
    }
    static int[] divide_pow_bin_(int[] a, int[] b) //V
    {
        int pow_a=0;
        int pow_b=0;
        int btoi_a=binaryToInt(a);
        int btoi_b=binaryToInt(b);
        ;
        while (power(2,pow_a)!=btoi_a)
            pow_a++;
        while (power(2,pow_b)!=btoi_b)
            pow_b++;
        int pow_res = Q-1+((pow_a-pow_b)%(Q-1));
        return power_bin_(new int[]{0,1},pow_res);
    }
    static int[] div_bin_(int[] a, int[] b)
    {return divide_bin_(a,b,new int[]{0})[0];}
    static int[] mod_bin_(int[] a, int[] b)
    {return divide_bin_(a,b,new int[]{0})[1];}

    static int[] temp_mult_bin(int[] a, int[] b) //V
    {
        int n = a.length - 1;
        int m = b.length - 1;
        int[] result = new int[n + m + 1];
        for (int i = 0; i <= n; i++) {
            for (int j = 0; j <= m; j++) {
                result[i + j] += a[i] * b[j];
                result[i + j] %= 2;
            }
        }
        return result;
    }
    static int[] multiply_bin_(int[] a, int[] b) //V
    {
        int n = a.length - 1;
        int m = b.length - 1;
        int[] result = new int[n + m + 1];
        for (int i = 0; i <= n; i++) {
            for (int j = 0; j <= m; j++) {
                result[i + j] += a[i] * b[j];
            }
        }

        for (int k=0; k<result.length; k++)
            result[k] %= 2;

        if (binaryToInt(result)<Q)
            return result;
        else return mod_bin_(result,intToBinary(Q_primal));
    }

    static int[] power_bin_(int[] a, int p)
    {
        int[] result = {1};
        for (int i = 0; i < p; i++)
            result=multiply_bin_(result, a);
        if (binaryToInt(result)<Q)
            return result;
        else return mod_bin_(result,intToBinary(Q_primal));
    }



    public static int plus_minus(int a, int b)
    {return binaryToInt(plus_minus_bin_(intToBinary(a),intToBinary(b)));}
    public static int divide(int a, int b)
    {return binaryToInt(div_bin_(intToBinary(a),intToBinary(b)));}
//    public static int modulo(int a, int b)
//    {return binaryToInt(mod_bin_(intToBinary(a),intToBinary(b)));}
    public static int multiply(int a, int b)
    {return binaryToInt(multiply_bin_(intToBinary(a),intToBinary(b)));}
    public static int power(int a, int p)
    {return binaryToInt(power_bin_(intToBinary(a),p));}
    public static int divpowbin(int a, int b)
    {return binaryToInt(divide_pow_bin_(intToBinary(a),intToBinary(b)));}

    public static int[] add_subPolynomialsGF(int[] a, int[] b) //V
    {
        int[] result = new int[Math.max(a.length, b.length)];
        for(int i=0; i<result.length; i++){
            int coeffA = (i<a.length) ? a[i] : 0;
            int coeffB = (i<b.length) ? b[i] : 0;
            result[i] = plus_minus(coeffA, coeffB);
        }
        return result;
    }

    public static int[] multiplyPolynomialsGF(int[] a, int[] b) //V
    {
        int[] result = new int[a.length + b.length - 1];
        for(int i=0; i<a.length; i++)
            for(int j=0; j<b.length; j++)
                result[i+j] = plus_minus(result[i+j], multiply(a[i], b[j]));
        result=Commons.trim(result);

        //if (binaryToInt(result)<Q)
            return result;
        //else
        //return modulePolynomials_GF(result,intToBinary(Q_primal));
    }
    public static int[][] dividePolynomialsGF_internal(int[] a, int[] b) //V
    {
        int m = a.length - 1;
        int n = b.length - 1;
        if (n == 0 || m < n)
            throw new IllegalArgumentException("Invalid input polynomial(s)");
        int[] q = new int[m - n + 1];
        int[] r = a;
        for (int i = m; i >= n; i--) {
            int coeff = divpowbin(r[i], b[n]);
            q[i - n] = coeff;
            for (int j = n; j >= 0; j--)
                r[i - n + j] = plus_minus(r[i - n + j],multiply(coeff, b[j]));
        }
        r = Commons.trim(r);
        return new int[][] {q, r};
    }
    public static int[] dividePolynomials_GF(int[] polyDividend, int[] polyDivisor)
    {return dividePolynomialsGF_internal(polyDividend, polyDivisor)[0];}
    public static int[] modulePolynomials_GF(int[] polyDividend, int[] polyDivisor)
    {return dividePolynomialsGF_internal(polyDividend, polyDivisor)[1];}

    //TODO
    public static int[] lagrangePolynomialsGF(int[] x, int[] y)
    {
        assert (x.length==y.length);

        int[] res = new int[]{0};
        for (int i = 0; i < x.length; i++)
            if (y[i] != 0) {
                int[] member = new int[]{1};
                int coeff = 1;
                for (int j = 0; j < x.length; j++)
                    if (i != j) {
                        int c = plus_minus(x[i], x[j]);
                        System.out.println(x[i] + " XOR " + x[j] + " = " + c);
                        coeff = multiply(coeff,c);
                    }
                System.out.println("coeff = " + coeff);


                for (int j = 0; j < x.length; j++)
                    if (i != j) {
                        int[] x_minus_xj = new int[]{Q-x[j],1};
                        //int[] x_minus_xj = new int[]{x[j],1};
                        System.out.println(Arrays.toString(member) + " *= " + Arrays.toString(x_minus_xj));
                        member = multiplyPolynomialsGF(member, x_minus_xj);
                        System.out.println(Arrays.toString(member));
                        //int coeff = plus_minus(x[i], x[j]);
                        //System.out.println(x[i] + " XOR " + x[j] + " = " + coeff);
                        //System.out.println(Arrays.toString(member) + " /= " + coeff);
//                        for (int k = 0; k < member.length; k++)
//                            //member[k] = divide(member[k], coeff);
//                            member[k] = divpowbin(member[k], coeff);
//                        System.out.println(Arrays.toString(member));
                    }

                System.out.println(Arrays.toString(member) + " /= " + coeff);
                for (int k = 0; k < member.length; k++)
                    member[k] = divpowbin(member[k], coeff);

                System.out.println(Arrays.toString(member) + " *= " + y[i]);
                for (int k = 0; k < member.length; k++)
                    member[k] = multiply(member[k], y[i]);

                System.out.println("Subtotal: " + Arrays.toString(member));

                System.out.println(Arrays.toString(res) + " XOR " + Arrays.toString(member));
                res = add_subPolynomialsGF(res, member);
                System.out.println("res subtotal: " + Arrays.toString(res));
            }
        return res;
    }
    public static int evalPolynomialGF(int x, int[] p) //V
    {
        int res=0;
        for (int i=0; i<p.length; i++)
        {
            int pw = power(x,i);
            System.out.println(x + " ** " + i + " = " + pw);

            int a = multiply(p[i],power(x,i));
            System.out.println(pw + " * " + p[i] + " = " + a);

            System.out.println(res + " XOR " + a + " = " + plus_minus(res, a));
            res=plus_minus(res,a);
            System.out.println();
        }
        if (res<Q)
            return res;
        else return binaryToInt(mod_bin_(intToBinary(res),intToBinary(Q_primal)));
    }
    private static void applyErrors(int[] p)
    {
        Random r = new Random();
        int j = Math.abs(r.nextInt() % p.length);
        int pj=p[j];
        while(pj==p[j])
            p[j] = Math.abs(r.nextInt() % Q);
    } //V

    public static int[] gaoEncode(int[] polynome, int[] message)
    {
        int[] c = new int[message.length];
        for (int i=0;i<message.length;i++)
            c[i]=evalPolynomialGF(message[i],polynome);
        System.out.println("Before error addition:\n"+Arrays.toString(c));
        applyErrors(c);
        return c;
    }

    //TODO
    static int[] partial_gcd_GF(int[] q0, int[] q1, int stop)
    {
        int[] r0 = q0;
        int[] r1=q1;

        int[] s0=Task1.intToBinary(1);
        int[] s1=Task1.intToBinary(0);

        int[] t0=Task1.intToBinary(0);
        int[] t1=Task1.intToBinary(1);

        while (r1.length-1>stop)
        {
            int[] q = Task1.dividePolynomials_GF(q0,q1);

            int[] ts=s0;
            int[] tt=t0;
            r0=r1;
            s0=s1;
            t0=t1;
            r1=Task1.modulePolynomials_GF(q0,q1);
            s1 = add_subPolynomialsGF(ts, multiplyPolynomialsGF(q, s0));
            t1 = add_subPolynomialsGF(tt, multiplyPolynomialsGF(q, t0));
        }
        return r1;
    }

    //TODO
    public static int[] gaoDecode(int[] a, int[] b)
    {
        int[] res = new int[]{};
        int[] g0 = new int[] {1};
        for (int i=0; i<a.length; i++)
            g0 = multiply_bin_(g0,new int[]{1, -1*a[i]}); //res*=x-a[i]

        int[] g1 = lagrangePolynomialsGF(a, b);

        return res;
    }
}

class Task2 extends Task1
{
    final static int BASE = 64;
    final static int l = 3;
    static int[] m = new int[l];
    static Map<Integer, Integer> s = new HashMap<Integer, Integer>();

    static void encode(String fingerprint)
    {
        Random r = new Random();
        for (int i=0; i<BASE; i++)
            if (fingerprint.indexOf(i)!=-1)
                s.put(i, evalPolynomialGF(i,m));
            else s.put(i, r.nextInt() % BASE);
    }
}


public class Main {
    public static void main(String[] args) {
//        System.out.println(Task1.plus_minus(6,3));
//        System.out.println(Task1.multiply(6,3));
//        System.out.println(Task1.divide(6,3)); может избавиться от обчычного деления напрочь?!!
//        System.out.println(Task1.divpowbin(6,3));
//        System.out.println(Task1.divpowbin(3,6));
//        System.out.println(Task1.modulo(5,3));
//        System.out.println(Task1.modulo(3,5));
//        System.out.println(Task1.power(2,3));
//        System.out.println(Task1.power(1,3));
//        System.out.println(Task1.power(5,2));


        int[] a = new int[]{1, 2, 3};
        int[] b = new int[]{3, 6, 7};
        int[] l = Task1.lagrangePolynomialsGF(a,b);




        System.out.println(Arrays.toString(l)); //[2,0,1]??
        for (int i=0; i<a.length; i++)
            System.out.printf("l("+ a[i] + ") = " + Task1.evalPolynomialGF(a[i],l) + "\n");

//        int[] aa = Task1.multiplyPolynomialsGF(new int[]{2,1}, new int[] {3,1});
//        int[] bb = Task1.multiplyPolynomialsGF(new int[]{2,1}, new int[] {1,1});
//        System.out.println(Task1.partial_gcd_GF(aa,bb,1));

    }
}
