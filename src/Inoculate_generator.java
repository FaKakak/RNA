import java.util.Random;

public class Inoculate_generator {
    public static final int INOCUNUM = 10;
    Environment environment;

    public Inoculate_generator(Environment environment) {
        this.environment = environment;
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
        Random random = new Random();
        int x = random.nextInt(Environment.SIDE);
        int y = random.nextInt(Environment.SIDE);
        RNA newRNA = new RNA(seq);
        environment.cell_head[h][y][x].addAfter(newRNA);
    }
}
