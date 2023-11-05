import java.util.LinkedList;
import java.util.List;
import java.util.Random;
import static java.lang.Math.pow;

public class Unit_case_generator {
    public static final double PSBP = 0.9;
    public static final double PBB = 0.000001;
    public static final double PLMC = 0.000002;
    public static final double PAT = 0.1;
    public static final double PLT = 0.9;
    public static final double PMOVR = 0.01;
    public static final double PMOV = (PMOVR/2);
    public static final double PMF = 0.0002;
    public static final int TNSS = 1;
    public static final double PMD = 0.01;
    public static final double PFP = 0.01;
    public static final double PMFS = 0.9;
    public double FDMOV(RNA p) {
        return pow(p.length1 + p.length2, 1 / 3.0);
    }

    Environment environment;
    int SIDE;
    int[][] raw_arr;
    RNA[][][] cell_head;
    Random random = new Random();
    int h;
    public Unit_case_generator(Environment environment) {
        this.environment = environment;
        SIDE = Environment.SIDE;
        raw_arr=environment.raw_arr;
        cell_head=environment.cell_head;
    }

    public void unit_case(int h){
        this.h = h;

        List<Integer> xy_init = new LinkedList<>();
        for (int i = 0; i < Environment.CELLNUM; i++) {
            xy_init.add(i);
        }
        int xyIndex,xy,x,y;

        for (int i = 0; i < Environment.CELLNUM; i++) {
            // 随机选x,y，相当于xy_choose(void)
            xyIndex = random.nextInt(xy_init.size());
            xy = xy_init.get(xyIndex);
            xy_init.remove(xyIndex);
            x = xy%SIDE;
            y = xy/SIDE;

            raw(y, x);

            for(RNA p = cell_head[h][y][x].next;p!=cell_head[h][y][x];p=p.next){
                switch(random.nextInt(6)){
                    case 0:
                        case0(p,y,x);
                        p = fresh_unit(p,y,x);
                        break;
                    case 1:
                        p=case1(p,y,x);
                        break;
                    case 2:
                        case2(p,y,x);
                        p = fresh_unit(p,y,x);
                        break;
                    case 3:
                        case3(p,y,x);
                        p = fresh_unit(p,y,x);
                        break;
                    case 4:
                        case4(p,y,x);
                        p = fresh_unit(p,y,x);
                        break;
                    case 5:
                        int[] rotate = case5(p,y,x);
                        p = fresh_unit(p,y+rotate[0],x+rotate[1]);
                        break;
                    default:
                }
            }
        }
    }
    private RNA fresh_unit(RNA p, int y, int x){
        // 好像成了
        // 它的功能就是，把当前RNA剥离这一层，搬到下一层，返回值为它的前一个RNA
        // 游离的RNA单独处理
        RNA p1 = p.prior;
        p.removeThis();
        cell_head[h^1][y][x].addAfter(p);
        return p1;
    }
    private boolean findSeq(char[] subseq,RNA rna){
        boolean hasSeq = false;
        if(rna.length1>= subseq.length){
            for (int i = 0; rna.length1- subseq.length>=i; i++) {
                int j;
                for (j = 0; j < subseq.length; j++) {
                    if(rna.information[0][i+j]!=subseq[j]) break;
                }
                if(j == subseq.length){
                    hasSeq = true;
                    break;
                }
            }
        }
        return hasSeq;
    }
    private void raw(int y, int x){
        int raw_bef = raw_arr[y][x];
        for(int i=0;i<raw_bef;i++){
            switch (random.nextInt(2)){
                case 0: // 原材料变成核苷酸
                    if(random.nextDouble()<PMF){
                        raw_arr[y][x]--;
                        RNA newRNA = new RNA();
                        switch (random.nextInt(4)+1){
                            case 1:
                                newRNA.information[0][0]='A';
                                break;
                            case 2:
                                newRNA.information[0][0]='C';
                                break;
                            case 3:
                                newRNA.information[0][0]='G';
                                break;
                            case 4:
                                newRNA.information[0][0]='U';
                                break;
                        }
                        newRNA.information[0][1]='0';
                        newRNA.information[1][0]='0';
                        newRNA.length1=1;
                        newRNA.length2=0;
                        newRNA.nick=0;
                        cell_head[h^1][y][x].addAfter(newRNA);
                    }
                    break;
                case 1: // 原材料游走
                    if(random.nextDouble()<PMOVR)
                    {
                    int direction = random.nextInt(4);   // Four possible directions
                    switch(direction)
                    {
                        case 0: // 左
                            if(x>0) // 最左的原材料会碰壁
                            {
                                raw_arr[y][x]--;
                                raw_arr[y][x-1]++;
                            }
                            break;

                        case 1: // 右
                            if(x<SIDE-1) // 最右的原材料会碰壁
                            {
                                raw_arr[y][x]--;
                                raw_arr[y][x+1]++;
                            }
                            break;

                        case 2: // 上
                            if(y>0) //最上的原材料会碰壁
                            {
                                raw_arr[y][x]--;
                                raw_arr[y-1][x]++;
                            }
                            break;

                        case 3: // 下
                            if(y<SIDE-1) // 最下的原材料会碰壁
                            {
                                raw_arr[y][x]--;
                                raw_arr[y+1][x]++;
                            }
                            break;

                        default:
                            System.err.println("raw moving error");
                            System.exit(1);
                    }
                }
                break;
                default:
                    System.err.println("raw case error");
                    System.exit(1);
            }
        }
    }

    private void case0(RNA p, int y, int x){
        //Chain ligation
        // 貌似done了
        for(RNA p3 = p.next;p3!=p;p3=p3.next){
            if (p3 == cell_head[h][y][x]) { p3 = cell_head[h][y][x].next; if (p3 == p)break; }

            if (p3.length2 == 0)
            {
                if(random.nextDouble()<PLMC/p3.length1){
                    if(p.length1+p3.length1>RNA.MAX_RNA_LENGTH){
                        continue;
                    }

                    for (int a = 0; a < p3.length1; a++) {
                        p.information[0][a+p.length1]=p3.information[0][a];
                    }
                    p.information[0][p.length1+p3.length1]='0';
                    p.length1=p.length1+p3.length1;
                    p3.removeThis();
                    break;
                }
            }
        }
    }

    private RNA case1(RNA p, int y, int x){
        if (p.length1 == 1)  // Decay of mononucleotide
        {
            // 这可能有问题
            if (p.length2 == 0 && random.nextDouble() < PMD)
            {
                // 单核苷酸降解
                raw_arr[y][x]++;
                p = fresh_unit(p,y,x);
                cell_head[h^1][y][x].next.removeThis();
                return p;
            }
            p = fresh_unit(p,y,x);
        }
        else                  //Degradation of chain
        {
            double f = PBB;
            for (int j = p.length1; j > 1; j--)
            {
                if (j <= p.length2)   // Falling into double chain region 这里还没检查其实
                {
                    int m,n,k;
                    if (p.nick == 0)
                    {
                        m = j - 1;
                        n = p.length2 - j + 1;
                        k = Math.min(m, n);
                        f = PBB * k * PBB * k;
                    }
                    else
                    {
                        if (j == p.nick + 1)
                        {
                            m = j - 1;
                            n = p.length2 - j + 1;
                            k = Math.min(m, n);
                            f = PBB * k;
                        }
                        else if (j > p.nick + 1)
                        {
                            m = j - p.nick - 1;
                            n = p.length2 - j + 1;
                            k = Math.min(m, n);
                            f = PBB * k * PBB * k;
                        }
                        else
                        {
                            m = j - 1;
                            n = p.nick - j + 1;
                            k = Math.min(m, n);
                            f = PBB * k * PBB * k;
                        }
                    }
                }

                if (random.nextDouble() < f)
                {
                    RNA p3 = new RNA();

                    for (int b = 0; b < p.length1 - j + 1; b++)
                        p3.information[0][b] = p.information[0][b + j - 1];
                    p3.length1 = p.length1 - j + 1;
                    p.length1 = j - 1;
                    p.information[0][p.length1] = '0';
                    p3.information[0][p3.length1] = '0';

                    if (p.length2 > j - 1)
                    {
                        for (int b = 0; b < p.length2 - j + 1; b++)
                            p3.information[1][b] = p.information[1][b + j - 1];
                        p3.length2 = p.length2 - j + 1;
                        p.length2 = j - 1;
                        p.information[1][p.length2] = '0';
                        p3.information[1][p3.length2] = '0';
                    }
                    else
                    {
                        p3.length2 = 0;
                        p3.information[1][p3.length2] = '0';
                    }

                    if (p.nick > j - 1) { p3.nick = p.nick - j + 1; p.nick = 0; }
                    else if (p.nick == j - 1) { p3.nick = 0; p.nick = 0; }
                    else p3.nick = 0;

                    fresh_unit(p3,y,x);
                    break;
                }
            }
            p = fresh_unit(p,y,x);
        }
        return p;
    }

    private void case2(RNA p, int y, int x){
        // 好像成了
        if (p.nick == 0)                      //Template-directed attraction of substrates
        {
            for (RNA p3 = p.next; p3 != p; p3 = p3.next)
            {
                if (p3 == cell_head[h][y][x]) { p3 = cell_head[h][y][x].next; if (p3 == p)break; }
                if (p3.length2 == 0) // 两条单链互补
                {
                    if (p3.length1 <= p.length1 - p.length2)
                    {
                        int flag= 1;
                        for (int b = 0; b < p3.length1; b++)
                        {
                            int combine = p3.information[0][p3.length1 - 1 - b] + p.information[0][p.length2 + b];
                            if (combine == 'A'+'U'||combine=='C'+'G') {flag = 0; }
                            else if (random.nextDouble()< PFP) {flag = 0; }
                            else { flag = 1; break; }
                        }
                        if (flag == 0)
                        {
                            double rtdaddphili = random.nextDouble();
                            if (rtdaddphili < PAT)
                            {
                                for (int a = 0; a < p3.length1; a++)
                                    p.information[1][p.length2 + a] = p3.information[0][p3.length1 - 1 - a];
                                p.information[1][p.length2 + p3.length1] = '0';
                                if (p.length2 != 0)p.nick = p.length2;
                                p.length2 = p.length2 + p3.length1;

                                p3.removeThis();
                                break;
                            }
                        }
                    }
                }
            }
        }
        else                     //Template-directed ligation
        {
            double rtdaddlig = random.nextDouble();
            if (rtdaddlig < PLT)
            {
                p.nick = 0;
            }
        }
    }

    private void case3(RNA p, int y, int x){
        // 应该成了
        if (p.length2 != 0)    // Separation of double chain
        {
            if (random.nextDouble()< pow(PSBP, p.length2 - p.nick))
            {
                RNA p3 = new RNA();
                for (int b = 0; b < p.length2 - p.nick; b++)
                    p3.information[0][b] = p.information[1][p.length2 - 1 - b];
                p.information[1][p.nick] = '0';

                p3.information[0][p.length2 - p.nick] = '0';
                p3.information[1][0] = '0';

                p3.length1 = p.length2 - p.nick;
                p.length2 = p.nick;
                p.nick = 0;
                p3.length2 = 0;
                p3.nick = 0;

                fresh_unit(p3,y,x);
            }
        }
    }

    private void case4(RNA p, int y, int x){
        // done了吧
        if (p.length2 == 0)
        {
            boolean flag = findSeq(RNA.NRSEQ,p);
            if (flag)    // nt-synthetase catalyses the synthesis of nt.
            {
                int nt_turn = TNSS;
                int raw_bef = raw_arr[y][x];
                for (int k = 0; k < raw_bef; k++)
                {
                    if (nt_turn <= 0)break;
                    nt_turn--;
                    if (random.nextDouble() < PMFS)
                    {
                        raw_arr[y][x]--;

                        RNA p3 = new RNA();
                        int randnt = random.nextInt(4)+ 1;
                        switch (randnt)
                        {
                            case 1:  p3.information[0][0] = 'A'; break;
                            case 2:  p3.information[0][0] = 'C'; break;
                            case 3:  p3.information[0][0] = 'G'; break;
                            case 4:  p3.information[0][0] = 'U'; break;
                            default:
                        }

                        p3.information[0][1] = '0';
                        p3.information[1][0] = '0';

                        p3.length1 = 1;
                        p3.length2 = 0;
                        p3.nick = 0;

                        fresh_unit(p3,y,x);
                    }
                }
            }
        }
    }

    private int[] case5(RNA p, int y, int x){
        int[] rotate = new int[]{0,0};
        if (random.nextDouble()* FDMOV(p) < PMOV)
        {
            int randcase1 = random.nextInt(4);   // Four possible directions
            switch (randcase1)
            {
                case 0:
                    if (x > 0)
                    {
                        rotate[1]=-1;
                    }
                    break;

                case 1:
                    if (x < SIDE - 1)
                    {
                        rotate[1]=1;
                    }
                    break;

                case 2:
                    if (y > 0)
                    {
                        rotate[0]=-1;
                    }
                    break;

                case 3:
                    if (y < SIDE - 1)
                    {
                        rotate[0]=1;
                    }
                    break;

                default:
            }
        }
        return rotate;
    }
}