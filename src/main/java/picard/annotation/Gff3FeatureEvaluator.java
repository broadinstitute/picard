package picard.annotation;

import org.apache.commons.lang3.tuple.Pair;
import org.obolibrary.obo2owl.OWLAPIObo2Owl;
import org.obolibrary.oboformat.parser.OBOFormatConstants;

import org.semanticweb.HermiT.ReasonerFactory;
import org.semanticweb.owlapi.apibinding.OWLManager;
import org.semanticweb.owlapi.model.*;
import org.semanticweb.owlapi.reasoner.OWLReasoner;
import org.semanticweb.owlapi.search.EntitySearcher;
import org.semanticweb.owlapi.vocab.OWLRDFVocabulary;
import picard.PicardException;

import java.net.URL;
import java.util.HashMap;
import java.util.Map;

public class Gff3FeatureEvaluator {
    final private OWLOntology ontology;
    final private OWLReasoner reasoner;

    final private HumanReadableNameVisitor visitor;

    final private Map<Pair<String, String>, Boolean> isSubClassOfCache = new HashMap<>();

    final private OWLDataFactory df = OWLManager.getOWLDataFactory();

    final static URL DEFAULT_ONTOLOGY_URL = Gff3FeatureEvaluator.class.getResource("so.obo");

    public Gff3FeatureEvaluator() {
        this(DEFAULT_ONTOLOGY_URL);
    }

    public Gff3FeatureEvaluator(final URL ontologyIRI) {
        final OWLOntologyManager manager = OWLManager.createOWLOntologyManager();
        try {
            ontology = manager.loadOntology(IRI.create(DEFAULT_ONTOLOGY_URL));
            final ReasonerFactory reasonerFactory = new ReasonerFactory();
            reasoner = reasonerFactory.createReasoner(ontology);

            visitor = new HumanReadableNameVisitor(ontology);

            //load human readable names map
            ontology.classesInSignature().forEach(cls -> cls.accept(visitor));
        } catch (final OWLOntologyCreationException ex) {
            throw new PicardException("error loading onotology from " + ontologyIRI);
        }
    }

    public boolean isASubClassOf(final String type, final String superClassType) {
        return isSubClassOfCache.computeIfAbsent(Pair.of(type, superClassType), p -> decideIfSubClassOf(p.getLeft(), p.getRight()));
    }

    private boolean decideIfSubClassOf(final String type, final String superClassType) {
        final OWLClassExpression cls = visitor.getOWLClass(type);
        final OWLClassExpression superCls = visitor.getOWLClass(superClassType);
        if (cls == null) {
            throw new PicardException("type " + type + " used in GFF3 file not found in ontology.");
        }
        if(superCls == null) {
            throw new PicardException("type " + superClassType + " used in GFF3 file not found in ontology.");
        }
        final OWLAxiom isSubClassOfAxiom = df.getOWLSubClassOfAxiom(cls, superCls);
        return reasoner.isEntailed(isSubClassOfAxiom);
    }

    private static class HumanReadableNameVisitor implements OWLClassExpressionVisitor {
        final private OWLAnnotationProperty idAnnotationProperty;
        final private OWLAnnotationProperty labelAnnotationProperty;
        final private Map<String, OWLClassExpression> humanReadableToOWLClassMap = new HashMap<>();
        final OWLDataFactory df = OWLManager.getOWLDataFactory();
        final OWLOntology ontology;

        HumanReadableNameVisitor(final OWLOntology ontology) {
            idAnnotationProperty = df.getOWLAnnotationProperty(OWLAPIObo2Owl.trTagToIRI(OBOFormatConstants.OboFormatTag.TAG_ID.getTag()));
            labelAnnotationProperty = df.getOWLAnnotationProperty(OWLRDFVocabulary.RDFS_LABEL.getIRI());
            this.ontology = ontology;

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
