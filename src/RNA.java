import java.lang.System;
import java.util.Arrays;

public class RNA {
    public static final int MAX_RNA_LENGTH = 100;
    public static final char[] NRSEQ = {'A','G','C','A','U','G','C','U'};
    public static final char[] CTR1SEQ = {'C','U','C','U','A','G','A','G'};
    public static final char[] CTR2SEQ = {'C','G','U','U','A','A','C','G'};
    public static final char[] CTR3SEQ = {'A','U','C','G','C','G','A','U'};
    public static final char[] INOCUSEQ = NRSEQ;
    public static final char[] INOCUSEQ1 = CTR1SEQ;
    public static final char[] INOCUSEQ2 = CTR2SEQ;
    public static final char[] INOCUSEQ3 = CTR3SEQ;
    char[][] information;
    int length1;
    int length2;
    int nick;

    RNA next;
    RNA prior;

    // 新RNA next指向自己，prior指向null
    public RNA() {
        information = new char[2][MAX_RNA_LENGTH];
        length1=0;
        length2=0;
        nick=0;
        next = this;
        prior = null;
        information[0][0] = '0';
        information[1][0] = '0';
    }
    public RNA(char[] bond1) {
        information = new char[2][MAX_RNA_LENGTH];
        for (int i = 0; i < bond1.length; i++) {
            information[0][i] = bond1[i];
        }
        length1 = bond1.length;
        length2 = 0;
        nick = 0;
        next = this;
        prior = null;
        information[0][length1] = '0';
        information[1][0] = '0';
    }

    public void addAfter(RNA newRNA){// 把newRNA添加到该RNA后边
        RNA oldNext = next;
        next = newRNA;
        newRNA.prior = this;
        newRNA.next = oldNext;
        if(oldNext.prior!=null) oldNext.prior = newRNA;
    }
    public void removeThis(){// 把当前RNA从链表中摘出来，变成自由的RNA
        // 游离的RNA不用动
        if(prior==null) {
            if(next == this) return ;
            else{
                System.err.println("cannot remove the head of linkedlist");
                System.exit(1);
            }
        }
        if(next.prior!=null) {
            // 链表头下边不止一个RNA
            next.prior = prior;
        }
        prior.next = next;
        next = this;
        prior = null;
    }
}
