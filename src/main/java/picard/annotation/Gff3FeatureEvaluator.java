package picard.annotation;

import org.apache.commons.lang3.tuple.Pair;
import org.obolibrary.obo2owl.OWLAPIObo2Owl;
import org.obolibrary.obo2owl.Obo2OWLConstants;
import org.obolibrary.oboformat.parser.OBOFormatConstants;

import org.semanticweb.HermiT.ReasonerFactory;
import org.semanticweb.owlapi.apibinding.OWLManager;
import org.semanticweb.owlapi.model.*;
import org.semanticweb.owlapi.reasoner.OWLReasoner;
import org.semanticweb.owlapi.search.EntitySearcher;
import org.semanticweb.owlapi.vocab.OWLRDFVocabulary;
import picard.PicardException;


import java.io.File;
import java.net.URL;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;

public class Gff3FeatureEvaluator {
    final private OWLOntology ontology;
    final private OWLReasoner reasoner;

    final private HumanReadableNameVisitor visitor = new HumanReadableNameVisitor();

    final private Map<Pair<String, String>, Boolean> isSubClassOfCache = new HashMap<>();
    final private Map<Pair<String, String>, Boolean> isPartOfCache = new HashMap<>();



    final private OWLDataFactory df = OWLManager.getOWLDataFactory();
    final OWLObjectPropertyExpression partOf;

    final static URL DEFAULT_ONTOLOGY_URL = Gff3FeatureEvaluator.class.getResource("so.owl");
    final public static String DEFAULT_CDS_IRI = "http://purl.obolibrary.org/obo/SO_0000316";
    final public static String DEFAULT_TRANSCRIPT_IRI = "http://purl.obolibrary.org/obo/SO_0000673";
    final public static String DEFAULT_GENE_IRI = "http://purl.obolibrary.org/obo/SO_0000704";
    final public static String DEFAULT_PARTOF_IRI = "http://purl.obolibrary.org/obo/so#part_of";
    final public static String DEFAULT_EXON_IRI = "http://purl.obolibrary.org/obo/SO_0000147";

    public Gff3FeatureEvaluator() {
        this(DEFAULT_ONTOLOGY_URL, DEFAULT_PARTOF_IRI);
    }

    public Gff3FeatureEvaluator(final URL ontologyIRI, final String partOfIRI) {
        final OWLOntologyManager manager = OWLManager.createOWLOntologyManager();
        try {
            //ontology = manager.loadOntology(IRI.create("file:///" + new File(ontologyIRI).getAbsolutePath()));
            ontology = manager.loadOntology(IRI.create(DEFAULT_ONTOLOGY_URL));
            final ReasonerFactory reasonerFactory = new ReasonerFactory();
            reasoner = reasonerFactory.createReasoner(ontology);

            partOf = df.getOWLObjectProperty(IRI.create(partOfIRI));

            //load human readable names map
            ontology.classesInSignature().forEach(cls -> cls.accept(visitor));
        } catch (final OWLOntologyCreationException ex) {
            throw new PicardException("error loading onotology from " + ontologyIRI);
        }
    }

    public boolean isASubClassOf(final String type, final String superClassIRI) {
        return isSubClassOfCache.computeIfAbsent(Pair.of(type, superClassIRI), p -> decideIfSubClassOf(p.getLeft(), p.getRight()));
    }

    private boolean decideIfSubClassOf(final String type, final String superClassIRI) {
        final OWLClassExpression cls = visitor.getOWLClass(type);
        if (cls == null) {
            throw new PicardException("type " + type + " used in GFF3 file not found in ontology.");
        }
        final OWLAxiom isSubClassOfAxiom = df.getOWLSubClassOfAxiom(cls, df.getOWLClass(IRI.create(superClassIRI)));
        return reasoner.isEntailed(isSubClassOfAxiom);
    }

    public boolean isPartOf(final String type, final String superClassIRI) {
        return isPartOfCache.computeIfAbsent(Pair.of(type, superClassIRI), p -> decideIfSubClassOf(p.getLeft(), p.getRight()));
    }

    private boolean decideIfPartOf(final String type, final String superClassIRI) {
        final OWLClassExpression cls = visitor.getOWLClass(type);
        if (cls == null) {
            throw new PicardException("type " + type + " used in GFF3 file not found in ontology.");
        }
        final OWLAxiom isPartOfAxiom = df.getOWLSubClassOfAxiom(cls, df.getOWLObjectSomeValuesFrom(partOf, df.getOWLClass(IRI.create(superClassIRI))));
        return reasoner.isEntailed(isPartOfAxiom);
    }

    class HumanReadableNameVisitor implements OWLClassExpressionVisitor {
        final private OWLAnnotationProperty idAnnotationProperty;
        final private OWLAnnotationProperty labelAnnotationProperty;
        final private Map<String, OWLClassExpression> humanReadableToOWLClassMap = new HashMap<>();
        final OWLDataFactory df = OWLManager.getOWLDataFactory();

        HumanReadableNameVisitor() {
            idAnnotationProperty = df.getOWLAnnotationProperty(OWLAPIObo2Owl.trTagToIRI(OBOFormatConstants.OboFormatTag.TAG_ID.getTag()));
            labelAnnotationProperty = df.getOWLAnnotationProperty(OWLRDFVocabulary.RDFS_LABEL.getIRI());

        }
        @Override
        public void visit(final OWLClass ce) {
            if (EntitySearcher.getAnnotationAssertionAxioms(ce, ontology).noneMatch(OWLAnnotationAssertionAxiom::isDeprecatedIRIAssertion)) {
                addIdentifierToMap(ce, idAnnotationProperty);
                addIdentifierToMap(ce, labelAnnotationProperty);
            }
        }

        void addIdentifierToMap(final OWLClass ce, final OWLAnnotationProperty annotationProperty) {
            EntitySearcher.getAnnotationObjects(ce, ontology, annotationProperty).forEach(an -> {
                        final String identifier = ((OWLLiteral)an.getValue()).getLiteral();
                        checkAndAddToMap(identifier, ce);
                    }
            );
        }

        void checkAndAddToMap(final String identifier, final OWLClass ce) {
            if (identifier != null) {
                if (humanReadableToOWLClassMap.containsKey(identifier)) {
                    throw new PicardException("Human readable identifier " + identifier + " found associated with multiple classes when loading GFF3 ontology.");
                }
                humanReadableToOWLClassMap.put(identifier, ce);
            }
        }

        OWLClassExpression getOWLClass(final String identifier) {
            return humanReadableToOWLClassMap.get(identifier);
        }
    }
}
