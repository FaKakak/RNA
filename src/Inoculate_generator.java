import java.util.Random;

public class Inoculate_generator {
    public static final int INOCUNUM = 10;
    Environment environment;
    Random random;
    public Inoculate_generator(Environment environment,Random random) {
        this.environment = environment;
        this.random = random;
    }

    public void inoculate(int h){
        for (int i = 0; i < INOCUNUM; i++) {
            inoculate_helper(h,RNA.INOCUSEQ);
            inoculate_helper(h,RNA.INOCUSEQ1);
            inoculate_helper(h,RNA.INOCUSEQ2);
            inoculate_helper(h,RNA.INOCUSEQ3);
        }
    }
    private void inoculate_helper(int h,char[] seq){
        int x = random.nextInt(Environment.SIDE);
        int y = random.nextInt(Environment.SIDE);
        RNA newRNA = new RNA(seq);
        environment.cell_head[h][y][x].addAfter(newRNA);
    }
}
