package picard.cmdline;

import javax.annotation.processing.AbstractProcessor;
import javax.annotation.processing.Filer;
import javax.annotation.processing.RoundEnvironment;
import javax.annotation.processing.SupportedAnnotationTypes;
import javax.lang.model.element.Element;
import javax.lang.model.element.TypeElement;
import javax.lang.model.type.MirroredTypeException;
import javax.lang.model.util.Elements;
import javax.tools.Diagnostic;
import javax.tools.FileObject;
import javax.tools.StandardLocation;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

/*
 * The MIT License
 *
 * Copyright (c) 2014 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
@SupportedAnnotationTypes("picard.cmdline.*")
public class CommandLineAnnotationProcessor extends AbstractProcessor {

    private class CommandLineNameComparator implements Comparator<String> {
        @Override
        public int compare(final String a, final String b) {
            return a.substring(a.lastIndexOf(".")+1).compareTo(b.substring(b.lastIndexOf(".")+1));
        }
    }

    @Override
    public boolean process(final Set<? extends TypeElement> annotations, final RoundEnvironment roundEnvironment) {
        if (roundEnvironment.processingOver()){
            return false;
        }

        final Map<String, Set<String>> services = new HashMap<String, Set<String>>();

        final Elements elements = processingEnv.getElementUtils();

        for (final Element element : roundEnvironment.getElementsAnnotatedWith(ProviderFor.class)) {
            final ProviderFor annotation = element.getAnnotation(ProviderFor.class);
            if (annotation == null) {
                continue;
            }
            if (!element.getKind().isClass() && !element.getKind().isInterface()) {
                continue;
            }
            String contractName = null;
            final TypeElement type = (TypeElement) element;
            try {
                annotation.value();
            } catch (final MirroredTypeException mte) {
                contractName = mte.getTypeMirror().toString();
            }

            if (contractName == null) {
                continue;
            }

            Set<String> allContractNames = services.get(contractName);
            if (allContractNames == null)
                services.put(contractName, allContractNames = new TreeSet<String>(new CommandLineNameComparator()));
            allContractNames.add(elements.getBinaryName(type).toString());
        }

        final Filer filer = processingEnv.getFiler();
        for (final Map.Entry<String, Set<String>> entry : services.entrySet()) {
            try {
                final String contract = entry.getKey();
                final FileObject f = filer.getResource(StandardLocation.CLASS_OUTPUT, "", "META-INF/services/" + contract);
                final BufferedReader r = new BufferedReader(new InputStreamReader(f.openInputStream(), "UTF-8"));
                String line;
                while ((line = r.readLine()) != null)
                    entry.getValue().add(line);
                r.close();
            } catch (final FileNotFoundException ignored) {
            } catch (final IOException x) {
                processingEnv.getMessager().printMessage(Diagnostic.Kind.ERROR, "Failed to load existing service definition files: " + x);
            }
        }

        for (final Map.Entry<String, Set<String>> entry : services.entrySet()) {
            try {
                final String contract = entry.getKey();
                processingEnv.getMessager().printMessage(Diagnostic.Kind.NOTE, "Writing META-INF/services/" + contract);
                final FileObject f = filer.createResource(StandardLocation.CLASS_OUTPUT, "", "META-INF/services/" + contract);
                final PrintWriter pw = new PrintWriter(new OutputStreamWriter(f.openOutputStream(), "UTF-8"));
                for (final String value : entry.getValue()) {
                    pw.println(value);
                }
                pw.close();
            } catch (final IOException x) {
                processingEnv.getMessager().printMessage(Diagnostic.Kind.ERROR, "Failed to write service definition files: " + x);
            } catch (final Exception e) {
                processingEnv.getMessager().printMessage(Diagnostic.Kind.ERROR, "Failed to write service definition files: " + e);
            }
        }
        return false;
    }
}
