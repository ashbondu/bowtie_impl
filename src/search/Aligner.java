package search;

import exceptions.NoChoiceAvailable;
import index.FMIndex;

import java.io.FileInputStream;
import java.io.ObjectInputStream;
import java.util.*;

/**
 * Created by ashwinsl on 12/1/15.
 */
public class Aligner {

    private static final int BACKTRACK_THRESHOLD = 1000;
    private static final int BACKTRACK_ITER_THRESHOLD = 50;
    private static final int HIGH_QUALITY_READ_VALUE = 100;
    private static final int READ_QUALITY_CUTOFF = 35;


    private Stack<Status> alignedTillNow = new Stack<>();
    private Status firstPositionStatus = null;

    //Using hashmap for now, can do a lot better.
    private HashMap<String,Object> map = new HashMap<>();
    private final String TOCHECK = "TOCHECK";
    private final String LOWINDEX = "LOWINDEX";
    private final String HIGHINDEX = "HIGHINDEX";

    private char[] dnaToAlign = null;
    private List<Integer> quality = null;

    private char[] lastCol;
    private int[][] tally;
    private int tallyLength;
    private HashMap<Integer, Integer> suffixArray = null;

    private int numberOfAs;
    private int numberOfCs;
    private int numberOfGs;
    private int numberOfTs;

    private int totalLength;


    public Aligner(String path){
        //Path is the absolute path where the index and LRow are present.
        //Can change this based on how the index is going to be.

        try {

            ObjectInputStream in = new ObjectInputStream(new FileInputStream(path + "/extradata"));
            int[] firstCol = (int[]) in.readObject();
            FMIndex.FMDISTANCE = in.read();
            int suffixDistance = (int) in.read();

            in.close();

            numberOfAs = firstCol[1] - firstCol[0];
            numberOfCs = firstCol[2] - numberOfAs - 1;
            numberOfGs = firstCol[3] - numberOfAs - numberOfCs - 1 ;
            numberOfTs = firstCol[4] - numberOfAs - numberOfCs - numberOfGs - 1;

            totalLength = numberOfAs + numberOfCs + numberOfGs + numberOfTs;

            tallyLength = totalLength/FMIndex.FMDISTANCE;

            in = new ObjectInputStream(new FileInputStream(path + "/lastColFile"));
            lastCol = (char[]) in.readObject();
            in.close();


            in = new ObjectInputStream(new FileInputStream(path + "/tally"));
            tally = (int[][]) in.readObject();
            in.close();

            in = new ObjectInputStream(new FileInputStream(path + "/sa"));
            suffixArray = (HashMap<Integer, Integer>) in.readObject();
            in.close();


        } catch (Exception e){
            e.printStackTrace();
        }
    }

    /**
     * This function is called when we train to align the dna read.
     * The Index of the genome is already loaded before this is called using hte constructor.
     * @param dnaToAlign The read of the DNA which needs to be aligned to the Genome.
     * @param qualityScores THe quality scores of the DNA read which needs to e aligned.
     * @return 1 : Successful Exit, -1 : Something went wrong !!
     */
    public int align(char[] dnaToAlign, List<Integer> qualityScores){

        this.dnaToAlign = dnaToAlign;
        this.quality = qualityScores;

        //This stack contains all the data regarding which position is aligned to which position.
        alignedTillNow.clear();

        int iteration = 0;
        int length = dnaToAlign.length;

        int lastRankLow = 0;
        int lastRankHigh = 0;

        ArrayList<Character> checkedTillNow = new ArrayList<>();

        char toCheck = dnaToAlign[length -1];
        char previous = ' ';

        Character[] firstPositionChar = {toCheck};

        firstPositionStatus = new Status(0,0,toCheck,Arrays.asList(firstPositionChar));

        boolean completed = true;
        int backtrackTimes = 0;

        if(length <= 1){
            //throw new Exception();
            System.out.println("Length is one or less, so cannot perform");
        }

        int lowIndex = 0, highIndex =0;

        try {
            for (iteration=0; iteration < length; iteration++) {

                //For position 0, the number of the same character we are trying to align forms the initial range.
                if(iteration == 0){
                    switch (toCheck) {
                        case 'a':
                        case 'A':
                            lowIndex = 1;
                            highIndex = numberOfAs;
                            break;
                        case 'c':
                        case 'C':
                            lowIndex = numberOfAs + 1;
                            highIndex = numberOfAs + numberOfCs;
                            break;
                        case 'g':
                        case 'G':
                            lowIndex = numberOfAs + numberOfCs + 1;
                            highIndex = numberOfAs + numberOfCs + numberOfGs;
                            break;
                        case 't':
                        case 'T':
                            lowIndex = numberOfAs + numberOfCs + numberOfGs + 1;
                            highIndex = numberOfAs + numberOfCs + numberOfGs + numberOfTs;
                            break;
                        default: throw new Exception("Not a valid alphabet");
                    }

                    previous = toCheck;
                    toCheck = dnaToAlign[length-1-iteration-1];
                    checkedTillNow.add(toCheck);
                    alignedTillNow.add(new Status(lowIndex,highIndex,previous,checkedTillNow));

                    checkedTillNow.clear();
                    continue;

                }

                int previousOffset = offsetForLRow(toCheck);

                //This will denote the lowest rank of the next character that the query can align to.
                lastRankLow = getRank(toCheck,lowIndex,true);
                int newLowIndex = previousOffset + lastRankLow + 1;

                //This denotes the highest rank of the next character that the query can align to.
                lastRankHigh = getRank(toCheck,highIndex,false);
                int newHighIndex = previousOffset + lastRankHigh;

                //If the range of the ranks of the next character is not greater than one means
                // there is no match. Hence we need to backtrack.
                if(lastRankLow >= lastRankHigh){
                    //This is where we need to backtrack.
                    backtrackTimes ++;

                    int numberOfStepsBackTracked = 0;
                    try {
                        //Backtracking.
                        numberOfStepsBackTracked = backTrack(iteration,toCheck);
                    } catch (NoChoiceAvailable e) {
                        //This is thrown only when we have backtracked for more than the
                        // threshold assigned to backtrack in the same position. Hence say no match.
                        System.out.println("No match found after backtracking more than " + BACKTRACK_ITER_THRESHOLD + " times from the same position.\n");
                        return -1;
                    } catch (Exception e){
                        //This is thrown when there is an unhandled exception.
                        //So report that there is a bug!!
                        System.out.println("No!! There is bug in the code.");
                        return -1;
                    }
                    iteration = iteration - numberOfStepsBackTracked;
                    toCheck = (char)map.get(TOCHECK);
                    lowIndex = (int)map.get(LOWINDEX);
                    highIndex = (int)map.get(HIGHINDEX);

                    if(backtrackTimes >= BACKTRACK_THRESHOLD){
                        //Tried backtracking more than 400 times. No match even after that.
                        System.out.println("No match Found even after " + BACKTRACK_THRESHOLD + " backtracks.\n");
                        return -1;
                    }

                } else {

                    if(length - 1 - (iteration + 1) < 0){
                        //We have reached the end of the DNA alignment.
                        //And Yay!! We have a match. So push the last letter and Break out of the loop.
                        alignedTillNow.add(new Status(lowIndex,highIndex,toCheck,null));
                        break;
                    }
                    previous = toCheck;
                    toCheck = dnaToAlign[length-1-iteration-1];
                    checkedTillNow.add(toCheck);
                    lowIndex = newLowIndex;
                    highIndex = newHighIndex;
                    alignedTillNow.add(new Status(lowIndex,highIndex,previous,checkedTillNow));
                    checkedTillNow.clear();
                }

            }
        } catch (Exception e) {
            completed = false;
            //e.printStackTrace();
        }

        //Now need to get the starting locations of the string which were aligned.
        //And pop off the stack to get all the elements which were aligned.

        if(completed) {
            StringBuilder stringBuilder = new StringBuilder();
            //stringBuilder.append(toCheck);
            int size = alignedTillNow.size();
            for (int iter = 0; iter < size; iter++) {
                Status status = alignedTillNow.pop();
                stringBuilder.append(status.toUse);
            }
            System.out.println("Aligned Against : " + stringBuilder);

            System.out.println("Suffixs are : ");

            //Find the suffixs for all the positions matched.
            for (int j = lowIndex; j <= highIndex; j++) {
                System.out.print(getSuffixPosition(j) + "\t");
            }
            System.out.println("\n");
        }
        return 1;
    }

    /**
     * Does a backtrack of the exact match till a point where we find a read which have the lowest quality.
     * Uses a different character at the lowest read quality position and backtracks the progress made until then to that position
     * @param iter : The current iteration the alignment is in.
     * @param toCheck : the character which is checked upon currently.
     * @return : the number of steps the backtracking has occurred.
     * @throws NoChoiceAvailable
     */
    private int backTrack(int iter, char toCheck) throws NoChoiceAvailable{

        int length = quality.size();

        //Need this when we are trying to get the second lowest, when the first lowest is not matching
        BackTrackStatus backTrackStatus = new BackTrackStatus();
        backTrackStatus.backtrackIter = 0;
        backTrackStatus.qualityReads = quality.toArray(new Integer[length]);

        //System.out.println("BACKTRACKING");

        boolean finished = false;
        Status lastStatus = null;

        while (!finished && backTrackStatus.backtrackIter < BACKTRACK_ITER_THRESHOLD) {
            backTrackStatus.backtrackIter++;

            //Get the positon where there is a read with lowest read quality.
            int position = quality.size() - getSmallestReadPosition(length -1 -iter, backTrackStatus.qualityReads) - 1;


            //If the same position is returned, then we need to add that the current to check is added to already tried.
            if (position == iter) {
                lastStatus = alignedTillNow.peek();
                //lastStatus.addToAlreadyTried(toCheck);
            } else if (position == 0) {

                //If the position returned is 0. then we need to start all over again.
                if(firstPositionStatus.getAlreadyTired().size() < 4){
                    toCheck = firstPositionStatus.getPreviouslyUncheckedCharsAndUpdate();
                    alignedTillNow.clear();
                } else {
                    //System.out.println("No!!!!");
                    backTrackStatus.qualityReads[length - 1] = HIGH_QUALITY_READ_VALUE;
                    continue;
                }

                //checkedTillNow.addAll(lastStatus.alreadyTried);
                map.put(TOCHECK, toCheck);
                map.put(LOWINDEX, lastStatus.lowIndex);
                map.put(HIGHINDEX, lastStatus.highIndex);

                return iter + 1;

            } else {
                lastStatus = alignedTillNow.get(position-1);
            }

            try {
                //Get a new Character to try for the next iteration.
                toCheck = lastStatus.getPreviouslyUncheckedCharsAndUpdate();
                int popIndex = iter - 1;

                int returnValue = 0;
                for (; popIndex >= position; popIndex--) {
                    alignedTillNow.pop();
                    returnValue++;
                }

                map.put(TOCHECK, toCheck);
                map.put(LOWINDEX, lastStatus.lowIndex);
                map.put(HIGHINDEX, lastStatus.highIndex);

                finished = true;
                return returnValue + 1;

            } catch (NoChoiceAvailable noChoiceAvailable) {
                //throw noChoiceAvailable;

                //Setting the qualityRead at a given position which have no more choices available as high
                //So the next time we check for minimum something else is returned.
                backTrackStatus.qualityReads[length - 1 - position] = HIGH_QUALITY_READ_VALUE;
            }
        }

        lastStatus = alignedTillNow.peek();

        try {
            toCheck = lastStatus.getPreviouslyUncheckedCharsAndUpdate();
            map.put(TOCHECK,toCheck);
            map.put(LOWINDEX,lastStatus.lowIndex);
            map.put(HIGHINDEX,lastStatus.highIndex);

        } catch (NoChoiceAvailable noChoiceAvailable) {
            backTrackStatus.qualityReads[length-iter] = HIGH_QUALITY_READ_VALUE;
        }

        return 1;

    }

    /**
     * Returns the position in the quality score array where the lowest quality read is present.
     * @param position The position after which we need to find the lowest quality read.
     * @param qualityReads The array of quality reads.
     * @return Returns the position where the lowest quality score occurs.
     */
    private int getSmallestReadPosition(int position, Integer[] qualityReads){
        //Allowing to backtrack upto 128 positions only.
        int length = qualityReads.length;
        if(length - position == 1){
            return position;
        }
        int offLength = Math.min(128,length - position - 1);

        int minPosition = position;
        int small = qualityReads[position];
        position = position + 1;

        for (int i = 0 ; i < offLength ; i++ ){
            if(small > qualityReads[position]){
                small = qualityReads[position];
                minPosition = position;

                if(small <= READ_QUALITY_CUTOFF){
                    return minPosition;
                }
            }
            position ++;
        }

        return minPosition;
    }


    /**
     * @param index The index of the Last column from where we need to find the suffix.
     * @return Returns the suffix positon for the corresponding Last column entry
     */
    private int getSuffixPosition(int index){
        Integer value = suffixArray.get(index);
        if(value != null){
            return value;
        }

        int iteration = 0;
        for ( ; value == null ;iteration++ ) {

            char lastColChar = lastCol[index];
            int rank = getRank(lastColChar,index,false);

            int checkOffset = offsetForLRow(lastColChar);

            index = rank + checkOffset;

            value = suffixArray.get(index);
        }

        value += iteration + 1;

        return value;

    }

    /**
     * @param toUse The character for which we need to find the starting offset.
     * @return Returns the starting offset for the character passed.
     */
    private int offsetForLRow(char toUse){
        switch (toUse) {
            case 'A':
                return 0;
            case 'C':
                return numberOfAs ;
            case 'G':
                return (numberOfAs+numberOfCs);
            case 'T':
                return (numberOfAs+numberOfCs+numberOfGs);
        }
        return -1;
    }

    /**Using the Tally of the FM Index, this returns the number of times the passed character is seen till the
     * index which is also passed. Basically the rank of the character at that particular tally index.
     * @param toCheck The character whose rank we need to know
     * @param index The index of the Last row where we want to know the rank of the character passed
     * @param isLow Boolean to suggest we are trying to find the rank for the lower bound of the range
     * @return Returns the rank of the character which is passed at the given index.
     */
    private int getRank(char toCheck, int index, boolean isLow){

        int position = 0;
        if(isLow){
            index = index -1;
        }
        position = (index)/FMIndex.FMDISTANCE;
        int stepsToWalk = index%FMIndex.FMDISTANCE;

        boolean isNextCheckPoint = false;

        if(index%FMIndex.FMDISTANCE > FMIndex.FMDISTANCE/2){
            if(! (index > tallyLength * FMIndex.FMDISTANCE)) {
                position = position + 1;
                stepsToWalk = FMIndex.FMDISTANCE - stepsToWalk;
                isNextCheckPoint = true;
            }
        }


        int checkPointCount = tally[FMIndex.getReleventRowNumber(toCheck)][position];

        int seen;
        if(isNextCheckPoint) {
            for (seen=0; stepsToWalk > 0; stepsToWalk--){
                if(toCheck == lastCol[index + stepsToWalk]){
                    seen++;
                }
            }
            return checkPointCount - seen;
        } else {
            for (seen=0; stepsToWalk > 0; stepsToWalk--){
                if(toCheck == lastCol[index - stepsToWalk + 1]){
                    seen++;
                }
            }
            return checkPointCount + seen;
        }

    }

    /**
     * This class is used as a Status store in the stack which maintains the alignment progress.
     */
    private class Status{
        private int lowIndex = 0;
        private int highIndex = 0;

        private List<Character> alreadyTried = null;
        private char toUse;

        public Status(int lowIndex, int highIndex, char toUse, List<Character> alreadyTried) {
            this.lowIndex = lowIndex;
            this.highIndex = highIndex;
            if(alreadyTried != null) {
                this.alreadyTried = new ArrayList<>(alreadyTried);
            }
            this.toUse = toUse;
        }


        public List<Character> getAlreadyTired(){
            return alreadyTried;
        }

        public char getPreviouslyUncheckedCharsAndUpdate() throws NoChoiceAvailable{

            if(alreadyTried.size() >= 4){
                throw new NoChoiceAvailable();
            }
            char[] dnaAlphabet = {'A','C','G','T'};
            int rand = new Random().nextInt(4);
            for(int i = 0; i < 4; i++){
                if(alreadyTried.contains(dnaAlphabet[rand])){
                    rand = (rand + 1) % 4 ;
                    continue;
                }
                else {
                    break;
                }
            }
            alreadyTried.add(dnaAlphabet[rand]);
            return dnaAlphabet[rand];
        }
    }

    private class BackTrackStatus{

        private Integer[] qualityReads;
        private int backtrackIter = 0;

    }


}
