package index;

/**
 * Created by ashwinsl on 12/1/15.
 */
public class FMIndex {

    public static int FMDISTANCE = 30;

    public static int getReleventRowNumber(char toGet){
        switch (toGet) {
            case 'A':
                return 0;
            case 'C':
                return 1;
            case 'G':
                return 2;
            case 'T':
                return 3;
        }
        return -1;
    }


}
