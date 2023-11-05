import java.util.Random;

public class Inoculate_generator {
    public static final int INOCUNUM = 1;
    Environment environment;
    int[][] raw_arr;
    RNA[][][] cell_head;

    public Inoculate_generator(Environment environment) {
        this.environment = environment;
        this.raw_arr = environment.raw_arr;
        this.cell_head = environment.cell_head;
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
        cell_head[h][y][x].addAfter(newRNA);
    }
}
