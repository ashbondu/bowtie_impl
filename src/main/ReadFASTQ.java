package main;

import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.features.FeatureInterface;
import org.biojava.nbio.core.sequence.features.QualityFeature;
import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.biojava.nbio.sequencing.io.fastq.*;
import search.Aligner;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

/**
 * Created by ashwinsl on 12/1/15.
 */
public class ReadFASTQ {

    public static void main(String[] args)
    {
        try {

            if(args.length != 2){
                System.out.println("");
            }

            FastqReader fastqReader = new SangerFastqReader();
            List<DNASequence> sequences = new LinkedList<DNASequence>();

            for (Fastq fastq : fastqReader.read(new File(args[1]))) {
                sequences.add(FastqTools.createDNASequenceWithQualityScores(fastq));
            }

            Aligner aligner = new Aligner(args[0]);
            List<Integer> quality = new ArrayList<>();

            //Reading all the sequences and their quality scores from the fastq file
            for(DNASequence sequence1 : sequences){

                List<FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound>> features = sequence1.getFeaturesByType("qualityScores");

                QualityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> qualityScores = (QualityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound>)features.get(0);

                for(Number num : qualityScores.getQualities()){
                    quality.add(num.intValue());
                }

                //Calling the alignment for all the reads in the fastq.
                long a = System.currentTimeMillis();
                aligner.align(sequence1.getSequenceAsString().toUpperCase().toCharArray(),quality);
                System.out.println("Time taken in milliseconds : " + (System.currentTimeMillis() - a));

                quality.clear();

            }

        }catch (Exception e){
            e.printStackTrace();
        }
    }
}
