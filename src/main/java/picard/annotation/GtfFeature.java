package picard.annotation;

import htsjdk.samtools.util.Interval;
import htsjdk.tribble.Feature;
import htsjdk.tribble.annotation.Strand;

import java.util.*;
import java.util.stream.Collectors;

public class GtfFeature implements Feature {

    private final String contig;
    private final String source;
    private final String type;
    private final int start;
    private final int end;
    private final Strand strand;
    private final int phase;
    private final Map<String, String> attributes;
    private final List<GtfFeature> parents = new ArrayList<>();
    private final List<GtfFeature> children = new ArrayList<>();

    public GtfFeature(final String contig, final String source, final String type, final int start, final int end, final Strand strand, final int phase, final Map<String, String> attributes) {
        this.contig = contig;
        this.source = source;
        this.type = type;
        this.start = start;
        this.end = end;
        this.phase = phase;
        this.strand = strand;
        this.attributes = attributes;
    }

    public String getSource() {
        return source;
    }

    @Override
    public int getEnd() {
        return end;
    }

    public Strand getStrand() {
        return strand;
    }

    public int getPhase() {
        return phase;
    }

    public String getType() {return type;}

    @Override
    public String getContig() {
        return contig;
    }

    @Override
    public int getStart() {
        return start;
    }

    public String getAttribute(final String key) {
        return attributes.get(key);
    }

    public Map<String, String> getAttributes() { return attributes;}

    public List<GtfFeature> getParents() {return parents;}

    public List<GtfFeature> getChildren() {return children;}

    public boolean hasParents() {return parents.size()>0;}

    public boolean hasChildren() {return children.size()>0;}

    public void addParent(final GtfFeature parent) {
        parents.add(parent);
    }

    public void addChild(final GtfFeature child) {
        children.add(child);
    }

    @Override
    public boolean equals(Object other) {
        if (!other.getClass().equals(GtfFeature.class)) {
            return false;
        }

        final GtfFeature otherGtfFeature = (GtfFeature) other;

        return (otherGtfFeature.getContig().equals(contig) &&
                otherGtfFeature.getSource().equals(source) &&
                otherGtfFeature.getType().equals(type) &&
                otherGtfFeature.getStart() == start &&
                otherGtfFeature.getEnd()== end &&
                otherGtfFeature.getStrand() == strand &&
                otherGtfFeature.getPhase() == phase &&
                otherGtfFeature.getAttributes().equals(attributes) &&
                otherGtfFeature.getParents().equals(parents) &&
                otherGtfFeature.getChildren().equals(children));
    }

    @Override
    public int hashCode() {
        return shallowHashCode() + children.stream().map(GtfFeature::shallowHashCode).reduce(0, Integer::sum) + parents.stream().map(GtfFeature::shallowHashCode).reduce(0, Integer::sum);
    }

    private int shallowHashCode() {
        return contig.hashCode() + source.hashCode() + type.hashCode() + new Integer(start).hashCode() + new Integer(end).hashCode() +
                strand.hashCode() + new Integer(phase).hashCode() + attributes.hashCode();
    }

    public Set<GtfFeature> flatten() {
        final HashSet<GtfFeature> features = new HashSet<>(Collections.singleton(this));

        for(final GtfFeature child : children) {
            features.addAll(child.flatten());
        }
        return features;
    }

    public Set<Interval> getRibosomalIntervals() {
        if (type == "rRNA") {
            return new HashSet<>(Arrays.asList(new Interval(contig, start, end)));
        } else {
            final Set<Interval> ribosomalIntervals = new HashSet<>();
            if (hasChildren()) {
                for (final GtfFeature child : children) {
                    ribosomalIntervals.addAll(child.getRibosomalIntervals());
                }
            }
            return ribosomalIntervals;
        }
    }
}
