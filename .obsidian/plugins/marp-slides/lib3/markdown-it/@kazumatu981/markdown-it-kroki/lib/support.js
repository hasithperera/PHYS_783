'use strict';

/**
 * Diagram Languages are supported by kroki.io
 */
const LANGUAGES = [
    'actdiag',
    'blockdiag',
    'bpmn',
    'bytefield',
    'c4plantuml',
    'dbml',
    'ditaa',
    'dot',
    'd2',
    'erd',
    'excalidraw',
    'graphviz',
    'mermaid',
    'nomnoml',
    'nwdiag',
    'packetdiag',
    'pikchr',
    'plantuml',
    'rackdiag',
    'seqdiag',
    'svgbob',
    'umlet',
    'vega',
    'vegalite',
    'wavedrom',
];

/**
 * Image formats are supported by kroki.io
 */
const IMG_FORMATS = [
    'png', 'svg', 'jpeg', 'pdf', 'base64'
];

module.exports = {
    lnaguages: LANGUAGES,
    imageFormats: IMG_FORMATS,
    /**
     * test whether `lang` is supported diagram language by kroki.io
     * @param {string} lang target language
     * @returns is supported
     */
    languageSupports: (lang) => LANGUAGES.includes(lang),
    /**
     * test whether `format` is supported image format by kroki.io
     * @param {string} format name of image format like 'png', 'svg', ... etc
     * @returns is supported
     */
    imageFormatSupports: (format) => IMG_FORMATS.includes(format)
};